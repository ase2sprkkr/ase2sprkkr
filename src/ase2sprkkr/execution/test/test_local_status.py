import pytest

from ase2sprkkr.execution import (
    LocalExecutor, RunRecord, RunState, Execution, WorkingDirectory,
    ExecutorEndpoint, FileManifest,
)


def test_disappeared_local_process_is_unknown(monkeypatch):
    executor = LocalExecutor()
    run = RunRecord(
        id="run", state=RunState.RUNNING, task="SCF",
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory("/tmp"), "12345"
        ),
        files=FileManifest(),
    )

    def missing_process(pid, signal):
        raise ProcessLookupError

    monkeypatch.setattr("ase2sprkkr.execution.local.os.kill", missing_process)

    (status,) = executor.status(run)

    assert status.state == RunState.UNKNOWN


def test_lost_is_not_a_run_state():
    with pytest.raises(ValueError):
        RunState("lost")
