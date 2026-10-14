from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    ExecutorPool,
    RunReference,
    RunState,
    RunStatus,
    RunStore,
)


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=0)
        self.exception = None

    def submit(self):
        self.reference = RunReference(executor="finishing-test", native_id="job-1")
        self.state = RunState.QUEUED

    def backend_finished(self):
        self.state = RunState.FINISHED

    def finalize(self):
        if self.callback:
            self.callback(self.result)

    def wait(self):
        self.backend_finished()
        self.finalize()
        return self.result

    def run(self):
        self.submit()
        return self.wait()

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass


class _Executor(Executor):
    name = "finishing-test"

    def __init__(self):
        self.runner = None
        self.status_value = None

    def create_run(self, reader, spec, callback=None):
        self.runner = _Runner(callback)
        return self.runner

    def status(self, *references):
        assert self.status_value is not None
        return tuple(
            RunStatus(
                reference=reference,
                state=self.status_value,
                return_code=0 if self.status_value == RunState.FINISHED else 1,
            )
            for reference in references
        )


def _run(tmp_path):
    executor = _Executor()
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"), executor_pool=ExecutorPool()
    )
    run = manager.create_run(
        executor,
        reader=object(),
        spec=SimpleNamespace(directory=str(tmp_path)),
        task="DOS",
        dataset="case",
    )
    run.submit()
    return manager, executor, run


def test_live_terminal_backend_state_is_finishing_until_callback(tmp_path):
    manager, executor, run = _run(tmp_path)
    executor.runner.backend_finished()

    assert run.state == RunState.FINISHING
    stored = manager.store.get(run.id)
    assert stored.state == RunState.FINISHING

    executor.runner.finalize()

    stored = manager.store.get(run.id)
    assert stored.state == RunState.FINISHED
    assert stored.backend_state == RunState.FINISHED


def test_refresh_keeps_backend_finished_run_active_until_finalization(tmp_path):
    manager, executor, run = _run(tmp_path)
    executor.status_value = RunState.FINISHED

    (refreshed,) = manager.refresh((run.record,))

    assert refreshed.state == RunState.FINISHING
    assert refreshed.backend_state == RunState.FINISHED
    assert refreshed.return_code == 0
    assert refreshed in manager.store.runs(active=True)


def test_refresh_remembers_failed_backend_outcome_while_finishing(tmp_path):
    manager, executor, run = _run(tmp_path)
    executor.status_value = RunState.FAILED

    (refreshed,) = manager.refresh((run.record,))

    assert refreshed.state == RunState.FINISHING
    assert refreshed.backend_state == RunState.FAILED
    assert refreshed.return_code == 1
