import io
from types import SimpleNamespace

from ase2sprkkr.execution.remote import _write_state, status_remote
from ase2sprkkr.tools.commands import remote as remote_command
from ase2sprkkr.tools.commands.remote import _finalize_failed_worker_start


def test_failed_worker_start_marks_prestart_run_failed(tmp_path):
    run_id = "run-1"
    _write_state(
        tmp_path,
        run_id,
        id=run_id,
        state="queued",
        directory=str(tmp_path),
        pid=None,
        created_at="2026-10-13T08:00:00+00:00",
        finished_at=None,
    )

    _finalize_failed_worker_start(
        str(tmp_path),
        {"id": run_id},
        OSError("cannot start worker"),
    )

    state = status_remote(tmp_path, run_id)
    assert state["state"] == "failed"
    assert "cannot start worker" in state["error"]
    assert state["finished_at"] is not None


def test_failed_worker_start_does_not_overwrite_started_run(tmp_path):
    run_id = "run-1"
    _write_state(
        tmp_path,
        run_id,
        id=run_id,
        state="running",
        directory=str(tmp_path),
        pid=1234,
        created_at="2026-10-13T08:00:00+00:00",
        finished_at=None,
    )

    _finalize_failed_worker_start(
        str(tmp_path),
        {"id": run_id},
        OSError("late error"),
    )

    state = status_remote(tmp_path, run_id)
    assert state["state"] == "running"
    assert state["pid"] == 1234
    assert state.get("error") is None


def test_resolve_mpi_command_preserves_integer_process_count(monkeypatch):
    monkeypatch.setattr(
        "ase2sprkkr.execution.mpi.find_default_mpi_launcher",
        lambda: ("mpirun",),
    )
    monkeypatch.setattr(remote_command.sys, "stdin", io.StringIO('{"mpi": 4}'))
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(SimpleNamespace(remote_action="resolve-mpi"), None)

    assert output == [{"mpi": ("mpirun", "-np", "4")}]
