import pytest

from ase2sprkkr.execution import LocalExecutor, RunReference, RunState


def test_disappeared_local_process_is_unknown(monkeypatch):
    executor = LocalExecutor()
    reference = RunReference(executor="local", native_id="12345")

    def missing_process(pid, signal):
        raise ProcessLookupError

    monkeypatch.setattr("ase2sprkkr.execution.local.os.kill", missing_process)

    (status,) = executor.status(reference)

    assert status.state == RunState.UNKNOWN


def test_lost_is_not_a_run_state():
    with pytest.raises(ValueError):
        RunState("lost")
