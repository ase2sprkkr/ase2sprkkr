from threading import Event
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    Executor,
    ExecutorPool,
    RunState,
    RunStatus,
    RunStore,
)

from ase2sprkkr.execution.manager import ExecutionManager


class _Runner:
    def __init__(self):
        self.state = RunState.NEW
        self.backend_id = None
        self.result = SimpleNamespace(return_code=None)
        self.exception = None
        self.failure_phase = None
        self.failure_message = None
        self._done = Event()

    def submit(self):
        self.state = RunState.QUEUED
        self.backend_id = "new-job"

    def wait(self):
        self._done.wait()
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED
        self._done.set()

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _RecoveringExecutor(Executor):
    name = "recover"

    def __init__(self, recovered):
        self.recovered = recovered
        self.recovery_calls = []

    def create_project_run(self, spec, *, run_id, callback=None):
        return _Runner()

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        self.recovery_calls.append((run_id, submission, submitted_at))
        if self.recovered is NotImplemented:
            raise NotImplementedError
        return self.recovered


def _run(manager, executor, tmp_path):
    return manager.create_run(
        executor,
        spec=SimpleNamespace(directory=str(tmp_path)),
        process_result=lambda outcome: outcome,
        task="DOS",
        dataset="case",
    )


def test_confirmed_absent_stale_submission_releases_namespace(tmp_path):
    executor = _RecoveringExecutor(None)
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        executor_pool=ExecutorPool(),
        submitting_grace=0,
    )
    abandoned = _run(manager, executor, tmp_path)
    manager.reserve(abandoned.record)
    replacement = _run(manager, executor, tmp_path)

    replacement.submit()

    old = manager.store.get(abandoned.id)
    assert old.state == RunState.FAILED
    assert "not found" in old.message
    assert replacement.record.state == RunState.QUEUED
    assert executor.recovery_calls[0][0] == abandoned.id
    assert executor.recovery_calls[0][2] is not None


@pytest.mark.parametrize(
    "recovered, grace, expected_state, expected_backend, recovery_calls",
    [
        pytest.param(
            RunStatus(state=RunState.QUEUED, backend_id="existing-job"),
            0, RunState.QUEUED, "existing-job", 1,
            id="recovered-backend-blocks-duplicate",
        ),
        pytest.param(
            NotImplemented, 0, RunState.SUBMITTING, None, 1,
            id="unsupported-recovery-keeps-namespace",
        ),
        pytest.param(
            None, 3600, RunState.SUBMITTING, None, 0,
            id="fresh-submitting-is-not-recovered",
        ),
    ],
)
def test_unresolved_or_recovered_submission_keeps_duplicate_protection(
    tmp_path, recovered, grace, expected_state, expected_backend, recovery_calls,
):
    executor = _RecoveringExecutor(recovered)
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        executor_pool=ExecutorPool(),
        submitting_grace=grace,
    )
    uncertain = _run(manager, executor, tmp_path)
    manager.reserve(uncertain.record)
    replacement = _run(manager, executor, tmp_path)

    with pytest.raises(DuplicateSubmissionError):
        replacement.submit()

    stored = manager.store.get(uncertain.id)
    assert stored.state == expected_state
    assert stored.execution.backend_id == expected_backend
    assert len(executor.recovery_calls) == recovery_calls
