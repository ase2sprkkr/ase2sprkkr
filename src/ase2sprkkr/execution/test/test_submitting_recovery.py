from threading import Event
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    ExecutionManager,
    Executor,
    ExecutorPool,
    RunState,
    RunStatus,
    RunStore,
)


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

    def create_persistent_run(self, spec, *, run_id, callback=None):
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


def test_recovered_submission_keeps_duplicate_protection(tmp_path):
    executor = _RecoveringExecutor(
        RunStatus(state=RunState.QUEUED, backend_id="existing-job")
    )
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        executor_pool=ExecutorPool(),
        submitting_grace=0,
    )
    uncertain = _run(manager, executor, tmp_path)
    manager.reserve(uncertain.record)
    replacement = _run(manager, executor, tmp_path)

    with pytest.raises(DuplicateSubmissionError):
        replacement.submit()

    recovered = manager.store.get(uncertain.id)
    assert recovered.state == RunState.QUEUED
    assert recovered.execution.backend_id == "existing-job"


def test_unsupported_recovery_never_releases_submission(tmp_path):
    executor = _RecoveringExecutor(NotImplemented)
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        executor_pool=ExecutorPool(),
        submitting_grace=0,
    )
    uncertain = _run(manager, executor, tmp_path)
    manager.reserve(uncertain.record)
    replacement = _run(manager, executor, tmp_path)

    with pytest.raises(DuplicateSubmissionError):
        replacement.submit()

    assert manager.store.get(uncertain.id).state == RunState.SUBMITTING


def test_fresh_submitting_record_is_not_recovered(tmp_path):
    executor = _RecoveringExecutor(None)
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        executor_pool=ExecutorPool(),
        submitting_grace=3600,
    )
    uncertain = _run(manager, executor, tmp_path)
    manager.reserve(uncertain.record)
    replacement = _run(manager, executor, tmp_path)

    with pytest.raises(DuplicateSubmissionError):
        replacement.submit()

    assert executor.recovery_calls == []
