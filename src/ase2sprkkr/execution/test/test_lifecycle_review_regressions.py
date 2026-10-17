from types import SimpleNamespace

import pytest

import ase2sprkkr.execution.manager as manager_module
from ase2sprkkr.execution import (
    ExecutionManager,
    FailurePhase,
    ResultFilesRecord,
    RunReference,
    RunState,
    RunStore,
)
from ase2sprkkr.execution.detached import _read_state, _write_state
from ase2sprkkr.execution.remote_psij import _ActiveRun, _RemotePsijAgent


def test_result_snapshot_sync_failure_is_terminal_finalization_failure(
    tmp_path, monkeypatch
):
    output = tmp_path / "result.out"
    output.write_text("result")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_files=(str(output),)),
    )
    store.save(record)
    record = manager.transition(record, state=RunState.FINALIZING)

    def fail_sync(_path):
        raise OSError("cannot synchronize result snapshot")

    monkeypatch.setattr(manager_module, "sync_file", fail_sync)

    with pytest.raises(OSError, match="cannot synchronize result snapshot"):
        manager.finalize(record, state=RunState.FINISHED, return_code=0)

    failed = store.get(record.id)
    assert failed.state == RunState.FINALIZATION_FAILED
    assert failed.failure_phase == FailurePhase.FINALIZATION
    assert failed.return_code == 0
    assert failed.finished_at is not None
    assert failed.result_files.captured_at is None
    assert "cannot synchronize result snapshot" in failed.message


def _psij_state(tmp_path, *, state=RunState.SUBMITTING, **changes):
    run_id = "run"
    run_dir = tmp_path / ".ase2sprkkr" / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    values = {
        "id": run_id,
        "state": state.value,
        "directory": str(tmp_path),
        "backend_reference": None,
        "submitted_at": "2026-10-17T08:00:00+00:00",
        "started_at": None,
        "finished_at": None,
        "cancel_requested": False,
    }
    values.update(changes)
    return _write_state(tmp_path, run_id, **values)


def test_psij_recovered_failed_job_exit_code_does_not_prove_execution_started(
    tmp_path, monkeypatch
):
    state = _psij_state(tmp_path)
    reference = RunReference(
        executor="psij:slurm",
        native_id="job-1",
        directory=str(tmp_path),
    )

    class Executor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            return SimpleNamespace(
                reference=reference,
                state=RunState.FAILED,
                return_code=1,
                message="scheduler failed before start",
                failure_phase=None,
            )

    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: Executor())

    assert agent._recover_submitting("run", state, {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.FAILED.value
    assert recovered["return_code"] == 1
    assert recovered["started_at"] is None
    assert recovered["failure_phase"] == FailurePhase.LAUNCH.value


class _FailedBeforeRunning:
    state = RunState.FAILED
    return_code = 1

    def wait(self):
        return SimpleNamespace(return_code=1)


def test_psij_terminal_exit_code_without_running_evidence_is_launch_failure(tmp_path):
    _psij_state(
        tmp_path,
        state=RunState.QUEUED,
        backend_reference={
            "executor": "psij:slurm",
            "native_id": "job-1",
            "computer": None,
            "directory": str(tmp_path),
            "input_file": None,
            "output_file": None,
            "error_file": None,
        },
    )
    agent = _RemotePsijAgent(tmp_path)
    active = _ActiveRun(_FailedBeforeRunning(), RunState.QUEUED)
    agent._active["run"] = active

    agent._finish("run", active, False)

    failed = _read_state(tmp_path, "run")
    assert failed["state"] == RunState.FAILED.value
    assert failed["return_code"] == 1
    assert failed["started_at"] is None
    assert failed["failure_phase"] == FailurePhase.LAUNCH.value
