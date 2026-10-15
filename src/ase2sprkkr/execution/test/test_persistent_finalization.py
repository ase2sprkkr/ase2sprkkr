from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    RunReference,
    RunState,
    RunStore,
)


class _FinalizationRunner:
    def __init__(self, callback):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
        self._exception = None

    @property
    def exception(self):
        return self._exception

    def submit(self):
        self.reference = RunReference(executor="test", native_id="job-1")
        self.state = RunState.QUEUED

    def fail_after_backend_completion(self):
        self.state = RunState.FINISHED
        self._exception = ValueError("result parsing failed")
        self.callback(self.result)

    def wait(self):
        if self._exception is not None:
            raise self._exception
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass


class _FinalizationExecutor(Executor):
    name = "test"

    def create_run(self, reader, spec, callback=None):
        self.runner = _FinalizationRunner(callback)
        return self.runner


def test_post_backend_failure_does_not_publish_finished(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FinalizationExecutor()
    run = manager.create_run(
        executor,
        reader=object(),
        spec=SimpleNamespace(directory=str(tmp_path)),
        task="DOS",
    )

    run.submit()
    executor.runner.fail_after_backend_completion()

    record = store.get(run.id)
    assert record.state == RunState.FAILED
    assert record.backend_state == RunState.FINISHED
    assert record.message == "result parsing failed"
