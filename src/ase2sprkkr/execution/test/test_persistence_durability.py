import os
import stat
import sys
from pathlib import Path

import pytest

import ase2sprkkr.execution._durability as durability
import ase2sprkkr.execution.detached as detached
import ase2sprkkr.execution.detached_local as detached_local
import ase2sprkkr.execution.manager as manager_module
from ase2sprkkr.execution.base import FailurePhase, RunState
from ase2sprkkr.execution.detached import run_detached_worker, status_detached, submit_detached
from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.manager import (
    ExecutionManager,
    ResultFilesRecord,
    RunStore,
    SubmissionKey,
)


pytestmark = pytest.mark.skipif(
    os.name != "posix",
    reason="strong directory fsync semantics are currently POSIX-only",
)


def test_durable_write_orders_file_fsync_before_replace_and_directory_fsync(
    tmp_path, monkeypatch
):
    path = tmp_path / "state.json"
    events = []
    real_fsync = durability.os.fsync
    real_replace = durability.os.replace

    def tracked_fsync(descriptor):
        kind = "directory" if stat.S_ISDIR(os.fstat(descriptor).st_mode) else "file"
        events.append(("fsync", kind))
        return real_fsync(descriptor)

    def tracked_replace(source, destination):
        events.append(("replace", Path(destination).name))
        return real_replace(source, destination)

    monkeypatch.setattr(durability.os, "fsync", tracked_fsync)
    monkeypatch.setattr(durability.os, "replace", tracked_replace)

    durability.durable_write_json(path, {"state": "running"})

    assert events == [
        ("fsync", "file"),
        ("replace", "state.json"),
        ("fsync", "directory"),
    ]


def test_durable_write_propagates_directory_barrier_failure(tmp_path, monkeypatch):
    path = tmp_path / "state.json"
    path.write_text("old")

    def fail_directory_fsync(_directory):
        raise OSError("directory fsync failed")

    monkeypatch.setattr(durability, "fsync_directory", fail_directory_fsync)

    with pytest.raises(OSError, match="directory fsync failed"):
        durability.durable_write_text(path, "new")

    # The rename may already have happened; the important contract is that the
    # caller was not told the transition was durably committed.
    assert path.read_text() == "new"
    assert not tuple(tmp_path.glob(".*.tmp"))


def test_durable_mkdir_syncs_each_new_directory_entry(tmp_path, monkeypatch):
    synced = []
    monkeypatch.setattr(
        durability,
        "fsync_directory",
        lambda directory: synced.append(Path(directory)),
    )

    durability.durable_mkdir(tmp_path / "one" / "two", parents=True)

    assert synced == [tmp_path, tmp_path / "one"]


def test_durable_mkdir_syncs_parent_when_directory_already_exists(
    tmp_path, monkeypatch
):
    directory = tmp_path / "created-by-lock-helper"
    directory.mkdir()
    synced = []
    monkeypatch.setattr(
        durability,
        "fsync_directory",
        lambda path: synced.append(Path(path)),
    )

    durability.durable_mkdir(directory, parents=True, exist_ok=True)

    assert synced == [tmp_path]


def test_run_store_uses_durable_json_writer(tmp_path, monkeypatch):
    calls = []
    real_write = manager_module.durable_write_json

    def tracked_write(path, data, **options):
        calls.append((Path(path), options))
        return real_write(path, data, **options)

    monkeypatch.setattr(manager_module, "durable_write_json", tracked_write)
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="SCF", directory=tmp_path)

    store.save(record)

    assert len(calls) == 1
    assert calls[0][0].name == f"{record.id}.json"
    assert calls[0][1]["mode"] == 0o600


def test_reserve_syncs_inputs_before_submitting_record(tmp_path, monkeypatch):
    input_file = tmp_path / "input.inp"
    input_file.write_text("input")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        submission=SubmissionKey(
            endpoint=ExecutorEndpoint(),
            directory=str(tmp_path),
            task="SCF",
        ),
        result_files=ResultFilesRecord(input_files=(str(input_file),)),
    )
    events = []

    monkeypatch.setattr(
        manager_module,
        "sync_file",
        lambda path: events.append(("sync", Path(path).name)),
    )
    real_transition = manager.transition

    def tracked_transition(current, **changes):
        if changes.get("state") == RunState.SUBMITTING:
            events.append(("state", "submitting"))
        return real_transition(current, **changes)

    monkeypatch.setattr(manager, "transition", tracked_transition)

    reserved = manager.reserve(record)

    assert reserved.state == RunState.SUBMITTING
    assert events[:2] == [("sync", "input.inp"), ("state", "submitting")]


def test_finalize_syncs_outputs_before_terminal_record(tmp_path, monkeypatch):
    output_file = tmp_path / "output.out"
    output_file.write_text("result")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_files=(str(output_file),)),
    )
    store.save(record)
    record = manager.transition(record, state=RunState.FINALIZING)
    events = []

    monkeypatch.setattr(
        manager_module,
        "sync_file",
        lambda path: events.append(("sync", Path(path).name)),
    )
    real_transition = manager.transition

    def tracked_transition(current, **changes):
        if changes.get("state") == RunState.FINISHED:
            events.append(("state", "finished"))
        return real_transition(current, **changes)

    monkeypatch.setattr(manager, "transition", tracked_transition)

    finalized = manager.finalize(record, state=RunState.FINISHED, return_code=0)

    assert finalized.state == RunState.FINISHED
    assert events[:2] == [("sync", "output.out"), ("state", "finished")]
    assert len(finalized.result_files.output_versions) == 1


def _submit_local_without_helper(tmp_path, monkeypatch, run_id, output):
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *_args: object())
    return submit_detached(
        {
            "id": run_id,
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [sys.executable, "-c", "print('done')"],
            "output_file": str(output),
        }
    )


def test_detached_submission_commits_input_then_request_then_state(tmp_path, monkeypatch):
    input_file = tmp_path / "input.inp"
    input_file.write_text("input")
    events = []
    real_sync = detached.sync_file
    real_write = detached.durable_write_json

    def tracked_sync(path):
        events.append("input")
        return real_sync(path)

    def tracked_write(path, data, **options):
        events.append(Path(path).name)
        return real_write(path, data, **options)

    def tracked_spawn(*_args):
        events.append("spawn")
        return object()

    monkeypatch.setattr(detached, "sync_file", tracked_sync)
    monkeypatch.setattr(detached, "durable_write_json", tracked_write)
    monkeypatch.setattr(detached, "_spawn_detached_worker", tracked_spawn)

    submit_detached(
        {
            "id": "ordered-submission",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": ["unused"],
            "input_file": str(input_file),
        }
    )

    assert events == ["input", "request.json", "state.json", "spawn"]


def test_detached_local_syncs_streams_before_terminal_state(tmp_path, monkeypatch):
    output = tmp_path / "calculation.out"
    run_id = "durable-finish"
    _submit_local_without_helper(tmp_path, monkeypatch, run_id, output)
    events = []
    real_sync = detached_local.sync_stream
    real_write = detached._write_state_unlocked

    def tracked_sync(stream):
        events.append("stream")
        return real_sync(stream)

    def tracked_write(workspace, current_run_id, **changes):
        if changes.get("state") in {"finished", "failed", "cancelled"}:
            events.append("terminal")
        return real_write(workspace, current_run_id, **changes)

    monkeypatch.setattr(detached_local, "sync_stream", tracked_sync)
    monkeypatch.setattr(detached, "_write_state_unlocked", tracked_write)

    run_detached_worker(tmp_path, run_id)

    assert "stream" in events
    assert events.index("stream") < events.index("terminal")
    assert status_detached(tmp_path, run_id)["state"] == "finished"


def test_detached_stream_sync_failure_is_finalization_failure(tmp_path, monkeypatch):
    output = tmp_path / "calculation.out"
    run_id = "stream-sync-failure"
    _submit_local_without_helper(tmp_path, monkeypatch, run_id, output)

    def fail_sync(_stream):
        raise OSError("stream fsync failed")

    monkeypatch.setattr(detached_local, "sync_stream", fail_sync)

    with pytest.raises(OSError, match="stream fsync failed"):
        run_detached_worker(tmp_path, run_id)

    state = status_detached(tmp_path, run_id)
    assert state["state"] == "failed"
    assert state["failure_phase"] == FailurePhase.FINALIZATION.value
    assert state["return_code"] == 0
