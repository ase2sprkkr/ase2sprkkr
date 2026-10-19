from ase2sprkkr.execution import (
    
    ExecutorEndpoint,
    RunState,
    RunStore,
)

from ase2sprkkr.execution.manager import ExecutionManager


def _record(manager, tmp_path):
    return manager.new_record(
        task="DOS",
        directory=tmp_path,
        endpoint=ExecutorEndpoint(),
    )


def test_stale_transition_cannot_regress_newer_lifecycle_state(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    original = _record(manager, tmp_path)
    submitting = manager.transition(original, state=RunState.SUBMITTING)
    running = manager.transition(
        submitting,
        state=RunState.RUNNING,
        execution=submitting.execution.submitted("123"),
    )

    merged = manager.transition(
        submitting,
        state=RunState.QUEUED,
    )

    assert merged.state == RunState.RUNNING
    assert merged.execution.backend_id == "123"
    assert manager.store.get(original.id) == merged
    assert running.execution.backend_id == "123"


def test_none_does_not_erase_known_completion_metadata(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    original = _record(manager, tmp_path)
    stale = manager.transition(original, state=RunState.RUNNING)
    finished = manager.transition(
        stale,
        state=RunState.FINISHED,
        return_code=0,
        message="complete",
    )

    merged = manager.transition(
        stale,
        state=RunState.RUNNING,
        return_code=None,
        message=None,
    )

    assert merged.state == RunState.FINISHED
    assert merged.return_code == 0
    assert merged.message == "complete"
    assert manager.store.get(original.id) == merged
    assert finished.state == RunState.FINISHED


def test_late_conflicting_final_state_does_not_replace_known_final_outcome(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    original = _record(manager, tmp_path)
    running = manager.transition(original, state=RunState.RUNNING)
    manager.transition(
        running,
        state=RunState.FINISHED,
        return_code=0,
    )

    merged = manager.transition(
        running,
        state=RunState.FAILED,
        return_code=1,
    )

    assert merged.state == RunState.FINISHED
    # A stale conflicting final update must not overwrite the precise outcome.
    assert merged.return_code == 0
