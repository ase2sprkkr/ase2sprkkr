from types import SimpleNamespace

import pytest

import ase2sprkkr.outputs.persistent as persistent_output
from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    RunReference,
    RunState,
    RunStore,
)
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    attach_persistent_run,
    load_persistent_result,
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

    def fail_after_backend_completion(self, *, backend_state=RunState.FINISHED):
        self.state = backend_state
        self.result.return_code = 0 if backend_state == RunState.FINISHED else 1
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


def _run(tmp_path):
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
    return store, manager, executor, run


def test_post_backend_failure_is_finalization_failure(tmp_path):
    store, _manager, executor, run = _run(tmp_path)

    executor.runner.fail_after_backend_completion()

    record = store.get(run.id)
    assert record.state == RunState.FINALIZATION_FAILED
    assert record.message == "result parsing failed"
    assert record.finished_at is not None
    assert record.result_files.captured_at is None
    assert not hasattr(record, "backend_state")


def test_parser_failure_after_failed_execution_is_still_finalization_failure(tmp_path):
    store, _manager, executor, run = _run(tmp_path)

    executor.runner.fail_after_backend_completion(
        backend_state=RunState.EXECUTION_FAILED
    )

    record = store.get(run.id)
    assert record.state == RunState.FINALIZATION_FAILED
    assert record.return_code == 1
    assert record.message == "result parsing failed"


def test_failed_reconstruction_terminates_finalization_and_can_be_retried(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="NO_SUCH_TASK", directory=tmp_path)
    record = manager.transition(record, state=RunState.FINALIZING)

    with pytest.raises(PersistentResultUnavailableError):
        attach_persistent_run(manager, record)

    failed = store.get(record.id)
    assert failed.state == RunState.FINALIZATION_FAILED
    assert failed.finished_at is not None
    assert "unknown persisted effective task" in failed.message

    retrying = manager.retry_finalization(failed)
    assert retrying.state == RunState.FINALIZING
    assert retrying.finished_at is None
    assert retrying.message is None

    with pytest.raises(PersistentResultUnavailableError):
        attach_persistent_run(manager, retrying)

    failed_again = store.get(record.id)
    assert failed_again.state == RunState.FINALIZATION_FAILED
    assert "unknown persisted effective task" in failed_again.message


def test_reader_reconstruction_programming_error_fails_finalization_and_reraises(
    tmp_path, monkeypatch
):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(record, state=RunState.FINALIZING)
    error = RuntimeError("reader constructor failed")

    class BrokenReader:
        def __init__(self, *_args, **_kwargs):
            raise error

    monkeypatch.setattr(
        persistent_output.reader_registry,
        "class_for_task",
        lambda _task_definition: BrokenReader,
    )

    with pytest.raises(RuntimeError, match="reader constructor failed") as exc_info:
        attach_persistent_run(manager, record)

    assert exc_info.value is error
    failed = store.get(record.id)
    assert failed.state == RunState.FINALIZATION_FAILED
    assert failed.message == "reader constructor failed"
    assert failed.finished_at is not None


def test_transient_attach_failure_keeps_finalization_retryable(tmp_path, monkeypatch):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(record, state=RunState.FINALIZING)

    def fail_attach(*_args, **_kwargs):
        raise OSError("temporary ssh failure")

    monkeypatch.setattr(manager, "attach", fail_attach)

    with pytest.raises(OSError, match="temporary ssh failure"):
        attach_persistent_run(manager, record, reader=object())

    current = store.get(record.id)
    assert current.state == RunState.FINALIZING
    assert current.finished_at is None


def test_failed_finalization_must_be_retried_before_result_load(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(record, state=RunState.FINALIZING)
    record = manager.fail_finalization(record, ValueError("download failed"))

    with pytest.raises(
        PersistentResultUnavailableError,
        match="retry finalization first",
    ):
        load_persistent_result(manager, record, reader=object())

    assert store.get(record.id).state == RunState.FINALIZATION_FAILED
