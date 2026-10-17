import asyncio
import json
from datetime import datetime, timedelta, timezone
from threading import Event, Thread

import pytest

import ase2sprkkr.execution.detached as detached
import ase2sprkkr.execution.detached_local as detached_local
from ..remote import (
    _run_directory,
    _write_state,
    cancel_remote,
    clean_remote_runs,
    list_remote_runs,
    run_remote_worker,
    status_remote,
    submit_remote,
)


def _request(workspace, run_id, output_name):
    return {
        "id": run_id,
        "executor": "local",
        "directory": str(workspace),
        "argv": ["true"],
        "input_file": None,
        "stdin": None,
        "output_file": output_name,
        "error_file": output_name + ".err",
    }


def test_remote_run_metadata_is_isolated_per_run(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *a: None)

    first = submit_remote(_request(workspace, "run-one", "one.out"))
    second = submit_remote(_request(workspace, "run-two", "two.out"))

    assert first["id"] == "run-one"
    assert second["id"] == "run-two"
    assert status_remote(workspace, "run-one")["output_file"].endswith("one.out")
    assert status_remote(workspace, "run-two")["output_file"].endswith("two.out")

    for run_id in ("run-one", "run-two"):
        run_dir = workspace / ".ase2sprkkr" / "runs" / run_id
        assert (run_dir / "request.json").is_file()
        assert (run_dir / "state.json").is_file()

    assert not (workspace / ".ase2sprkkr-run.json").exists()
    assert not (workspace / ".ase2sprkkr-request.json").exists()


def test_remote_run_listing_keeps_history(tmp_path):
    workspace = tmp_path / "workspace"
    for run_id, state in (("old", "finished"), ("new", "failed")):
        run_dir = _run_directory(workspace, run_id)
        run_dir.mkdir(parents=True)
        (run_dir / "request.json").write_text("{}")
        _write_state(
            workspace,
            run_id,
            id=run_id,
            state=state,
            created_at=(
                datetime(2026, 1, 1, tzinfo=timezone.utc)
                if run_id == "old"
                else datetime(2026, 1, 2, tzinfo=timezone.utc)
            ).isoformat(),
            finished_at=datetime(2026, 1, 3, tzinfo=timezone.utc).isoformat(),
        )

    assert [run["id"] for run in list_remote_runs(workspace)] == ["old", "new"]


def test_remote_clean_removes_only_completed_metadata(tmp_path):
    workspace = tmp_path / "workspace"
    artifact = workspace / "result.dat"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("keep me")

    now = datetime.now(timezone.utc)
    states = {
        "finished": ("finished", now - timedelta(days=40)),
        "recent": ("failed", now - timedelta(days=1)),
        "active": ("running", None),
    }
    for run_id, (state, finished_at) in states.items():
        run_dir = _run_directory(workspace, run_id)
        run_dir.mkdir(parents=True)
        (run_dir / "request.json").write_text(json.dumps({"id": run_id}))
        _write_state(
            workspace,
            run_id,
            id=run_id,
            state=state,
            created_at=(now - timedelta(days=50)).isoformat(),
            finished_at=None if finished_at is None else finished_at.isoformat(),
        )

    cleaned = clean_remote_runs(workspace, older_than=30 * 86400)

    assert cleaned == {"removed": ["finished"], "active": ["active"]}
    assert not _run_directory(workspace, "finished").exists()
    assert _run_directory(workspace, "recent").exists()
    assert _run_directory(workspace, "active").exists()
    assert artifact.read_text() == "keep me"


def test_remote_clean_without_age_removes_all_completed_runs(tmp_path):
    workspace = tmp_path / "workspace"
    for run_id, state in (("done", "cancelled"), ("active", "queued")):
        run_dir = _run_directory(workspace, run_id)
        run_dir.mkdir(parents=True)
        _write_state(
            workspace,
            run_id,
            id=run_id,
            state=state,
            created_at=datetime.now(timezone.utc).isoformat(),
            finished_at=(
                datetime.now(timezone.utc).isoformat()
                if state == "cancelled"
                else None
            ),
        )

    cleaned = clean_remote_runs(workspace)

    assert cleaned == {"removed": ["done"], "active": ["active"]}
    assert not _run_directory(workspace, "done").exists()
    assert _run_directory(workspace, "active").exists()


def test_remote_submit_rejects_duplicate_run_id(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *a: None)
    request = _request(workspace, "same", "out")

    submit_remote(request)
    with pytest.raises(RuntimeError, match="already exists"):
        submit_remote(request)


def test_remote_cancel_waits_for_launch_handoff(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    run_id = "cancel-race"
    request = _request(workspace, run_id, "out")
    run_dir = _run_directory(workspace, run_id)
    run_dir.mkdir(parents=True)
    (run_dir / "request.json").write_text(json.dumps(request))
    _write_state(
        workspace,
        run_id,
        id=run_id,
        state="submitting",
        directory=str(workspace),
        pid=None,
        finished_at=None,
    )

    process_start_entered = Event()
    allow_process_start = Event()
    terminated = Event()
    cancel_done = Event()
    worker_errors = []
    killed = []

    class FakeProcess:
        pid = 4321
        returncode = None

        async def wait(self):
            assert await asyncio.to_thread(terminated.wait, 2)
            self.returncode = -15
            return self.returncode

    async def fake_start_local_process(*args, **kwargs):
        process_start_entered.set()
        assert await asyncio.to_thread(allow_process_start.wait, 2)
        return FakeProcess()

    def fake_killpg(pid, sig):
        killed.append((pid, sig))
        if sig == 0:
            if terminated.is_set():
                raise ProcessLookupError
            return
        terminated.set()

    monkeypatch.setattr(detached_local, "start_local_process", fake_start_local_process)
    monkeypatch.setattr(detached.os, "killpg", fake_killpg)

    def worker():
        try:
            run_remote_worker(workspace, run_id)
        except BaseException as exc:
            worker_errors.append(exc)

    worker_thread = Thread(target=worker)
    worker_thread.start()
    assert process_start_entered.wait(2)

    cancel_result = []

    def cancel():
        cancel_result.append(cancel_remote(workspace, run_id))
        cancel_done.set()

    cancel_thread = Thread(target=cancel)
    cancel_thread.start()

    assert not cancel_done.wait(0.1)
    allow_process_start.set()

    cancel_thread.join(2)
    worker_thread.join(2)

    assert not cancel_thread.is_alive()
    assert not worker_thread.is_alive()
    assert worker_errors == []
    assert cancel_result[0]["state"] == "cancelled"
    assert status_remote(workspace, run_id)["state"] == "cancelled"
    assert killed and killed[0][0] == FakeProcess.pid
