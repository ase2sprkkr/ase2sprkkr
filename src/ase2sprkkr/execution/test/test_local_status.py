import pytest

from ase2sprkkr.execution import (
    LocalExecutor, RunRecord, RunState, Execution, WorkingDirectory,
    ExecutorEndpoint, FileManifest,
)


def test_disappeared_local_process_is_unknown(monkeypatch):
    from ase2sprkkr.execution import remote_service

    executor = LocalExecutor()
    handle = remote_service._encode_handle({
        "kind": "process",
        "pid": 12345,
        "start_time": "0",
        "receipt_id": "a" * 32,
        "socket_path": "/missing/socket",
    })
    run = RunRecord(
        id="run", state=RunState.RUNNING, task="SCF",
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory("/tmp"), handle
        ),
        files=FileManifest(),
    )

    monkeypatch.setattr(remote_service, "_process_outcome", lambda receipt_id: None)
    monkeypatch.setattr(remote_service, "_process_info", lambda pid: None)

    (status,) = executor.status(run)

    assert status.state == RunState.UNKNOWN


def test_lost_is_not_a_run_state():
    with pytest.raises(ValueError):
        RunState("lost")
