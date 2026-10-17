from threading import Event, Thread
from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionManager,
    ExecutionOutcome,
    Executor,
    ExecutorPool,
    FailurePhase,
    PostprocessingRunner,
    RunState,
    RunStatus,
    RunStore,
)


class _Executor(Executor):
    name = "postprocessing-test"

    def __init__(self, status):
        self.status_value = status

    def status(self, *runs):
        return tuple(
            RunStatus(
                backend_id=run.execution.backend_id,
                state=self.status_value,
                return_code=0 if self.status_value == RunState.FINISHED else 1,
                failure_phase=(
                    FailurePhase.EXECUTION
                    if self.status_value == RunState.FAILED
                    else None
                ),
            )
            for run in runs
        )


class _RawRunner:
    def __init__(self):
        self.result = ExecutionOutcome(return_code=0)
        self.state = RunState.NEW
        self.backend_id = "job-1"
        self.exception = None
        self.failure_phase = None
        self._done = Event()

    def submit(self):
        self.state = RunState.RUNNING

    def wait(self):
        self._done.wait()
        self.state = RunState.FINISHED
        return self.result

    def complete(self):
        self._done.set()

    def cancel(self):
        self.state = RunState.CANCELLED
        self._done.set()

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


def _record(tmp_path, executor, state=RunState.RUNNING):
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"), executor_pool=ExecutorPool()
    )
    endpoint = manager.register_executor(executor)
    record = manager.new_record(
        task="DOS",
        directory=str(tmp_path),
        endpoint=endpoint,
    )
    record = manager.transition(
        record,
        state=state,
        execution=record.execution.submitted("job-1"),
    )
    return manager, record


def test_refresh_maps_raw_success_to_postprocessing(tmp_path):
    executor = _Executor(RunState.FINISHED)
    manager, record = _record(tmp_path, executor)

    (refreshed,) = manager.refresh((record,))

    assert refreshed.state == RunState.POSTPROCESSING
    assert refreshed.return_code == 0
    assert refreshed in manager.store.runs(active=True)


def test_refresh_keeps_raw_execution_failure_terminal(tmp_path):
    executor = _Executor(RunState.FAILED)
    manager, record = _record(tmp_path, executor)

    (refreshed,) = manager.refresh((record,))

    assert refreshed.state == RunState.FAILED
    assert refreshed.failure_phase == FailurePhase.EXECUTION
    assert refreshed.return_code == 1
    assert refreshed in manager.store.runs(active=False)


def test_postprocessing_is_observable_until_processor_finishes():
    raw = _RawRunner()
    entered = Event()
    release = Event()

    def process(outcome):
        entered.set()
        release.wait()
        return SimpleNamespace(return_code=outcome.return_code)

    runner = PostprocessingRunner(raw, process)
    runner.submit()
    raw.complete()

    assert entered.wait(1)
    assert runner.state == RunState.POSTPROCESSING

    release.set()
    result = runner.wait()

    assert result.return_code == 0
    assert runner.state == RunState.FINISHED
