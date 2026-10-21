from threading import Event
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    ExecutionOutcome,
    Executor,
    ExecutorPool,
    FailurePhase,
    RunState,
    RunStatus,
    RunStore,
)

from ase2sprkkr.execution.manager import ExecutionManager


class _RawRunner:
    def __init__(self, on_submit=None):
        self.on_submit = on_submit
        self.state = RunState.NEW
        self.backend_id = None
        self.result = ExecutionOutcome(return_code=None)
        self.exception = None
        self.failure_phase = None
        self.failure_message = None
        self._done = Event()

    @property
    def return_code(self):
        return self.result.return_code

    def submit(self):
        if self.on_submit:
            self.on_submit()
        self.backend_id = "native-1"
        self.state = RunState.QUEUED

    def wait(self):
        self._done.wait()
        return self.result

    def finish(self, return_code=0):
        self.result.return_code = return_code
        if return_code == 0:
            self.state = RunState.FINISHED
            self.failure_phase = None
            self.failure_message = None
        else:
            self.state = RunState.FAILED
            self.failure_phase = FailurePhase.EXECUTION
            self.failure_message = f"process exited with {return_code}"
        self._done.set()

    def cancel(self):
        self.state = RunState.CANCELLED
        self._done.set()

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _Executor(Executor):
    name = "test"

    def __init__(self):
        self.runners = []
        self.status_calls = []
        self.on_submit = None

    def create_project_run(self, spec, *, run_id, callback=None):
        runner = _RawRunner(on_submit=self.on_submit)
        self.runners.append(runner)
        return runner

    def status(self, *runs):
        self.status_calls.append(runs)
        return tuple(RunStatus(RunState.RUNNING, backend_id=run.execution.backend_id) for run in runs)


def _spec(tmp_path, *, dataset="case"):
    return SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(tmp_path / f"{dataset}.inp"),
        output_file=str(tmp_path / f"{dataset}.out"),
        input_files=(),
        output_files=(),
    )


def _manager(tmp_path):
    return ExecutionManager(
        RunStore(tmp_path / "runs"), executor_pool=ExecutorPool()
    )


def test_submit_persists_submitting_before_backend_call(tmp_path):
    manager = _manager(tmp_path)
    executor = _Executor()
    run = manager.create_run(
        executor,
        _spec(tmp_path),
        process_result=lambda outcome: outcome,
        task="ARPES",
        dataset="case",
    )
    seen = []
    executor.on_submit = lambda: seen.append(manager.store.get(run.id).state)
    # The runner already captured the callable during construction.
    executor.runners[0].on_submit = executor.on_submit

    run.submit()

    assert seen == [RunState.SUBMITTING]
    stored = manager.store.get(run.id)
    assert stored.state == RunState.QUEUED
    assert stored.finished_at is None
    assert stored.execution.backend_id == executor.runners[0].backend_id
    assert stored.task == "ARPES"
    assert stored.dataset == "case"


@pytest.mark.parametrize(
    "return_code, expected_state, expected_phase",
    [
        pytest.param(0, RunState.FINISHED, None, id="postprocess-success"),
        pytest.param(7, RunState.FAILED, FailurePhase.EXECUTION,
                     id="skip-processing-on-failure"),
    ],
)
def test_completion_postprocesses_only_successful_execution(
    tmp_path, return_code, expected_state, expected_phase,
):
    manager = _manager(tmp_path)
    executor = _Executor()
    processed = []
    run = manager.create_run(
        executor,
        _spec(tmp_path),
        process_result=lambda outcome: processed.append(outcome) or outcome,
        task="DOS",
    )
    run.submit()
    executor.runners[0].finish(return_code)

    result = run.wait()

    assert result.return_code == return_code
    assert processed == ([result] if return_code == 0 else [])
    stored = manager.store.get(run.id)
    assert stored.state == expected_state
    assert stored.failure_phase == expected_phase
    if return_code == 0:
        assert stored.finished_at is not None
    else:
        assert stored.return_code == return_code


def test_refresh_batches_status_queries(tmp_path):
    manager = _manager(tmp_path)
    executor = _Executor()
    endpoint = manager.register_executor(executor)
    records = []
    for index in range(2):
        record = manager.new_record(
            task="DOS",
            directory=tmp_path / str(index),
            endpoint=endpoint,
        )
        record = manager.transition(
            record,
            state=RunState.QUEUED,
            execution=record.execution.submitted(str(index)),
        )
        records.append(record)

    refreshed = manager.refresh(records)

    assert [record.state for record in refreshed] == [
        RunState.RUNNING,
        RunState.RUNNING,
    ]
    assert executor.status_calls == [tuple(records)]


def test_duplicate_submission_is_blocked_by_persisted_namespace(tmp_path):
    manager = _manager(tmp_path)
    executor = _Executor()
    first = manager.create_run(
        executor,
        _spec(tmp_path),
        process_result=lambda outcome: outcome,
        task="DOS",
        dataset="case",
    )
    first.submit()

    second = manager.create_run(
        executor,
        _spec(tmp_path),
        process_result=lambda outcome: outcome,
        task="DOS",
        dataset="case",
    )

    with pytest.raises(DuplicateSubmissionError):
        second.submit()
