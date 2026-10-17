from types import SimpleNamespace

import pytest

from ase2sprkkr.common.process_output_reader import OutputParseError
from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    ResultFilesRecord,
    RunReference,
    RunState,
    RunStore,
)
from ase2sprkkr.outputs.persistent import (
    PersistentResultMismatchError,
    PersistentResultUnavailableError,
    UnrecoverablePersistentResultError,
    check_persistent_result,
    load_persistent_result,
)


class _ResultReader:
    print_output = False

    def __init__(self, store=None, *, errors=()):
        self.store = store
        self.run_id = None
        self.errors = list(errors)
        self.calls = 0

    def create_result(self, output_file, input_file=None):
        return SimpleNamespace(
            output_file=output_file,
            input_file=input_file,
            return_code=None,
            complete=lambda *_args: None,
        )

    def read_from_file(self, filename, **_kwargs):
        self.calls += 1
        if self.store is not None and self.run_id is not None:
            assert self.store.get(self.run_id).state == RunState.AVAILABLE
        if self.errors:
            raise self.errors.pop(0)
        return SimpleNamespace(filename=filename, run_record=None)


class _LifecycleRunner:
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

    def complete(self, state=RunState.FINISHED, *, error=None):
        self.state = state
        self.result.return_code = 0 if state == RunState.FINISHED else 1
        self._exception = error
        self.callback(self.result)

    def wait(self):
        if self._exception is not None:
            raise self._exception
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _LifecycleExecutor(Executor):
    name = "test"

    def create_run(self, reader, spec, callback=None):
        self.capture_reader = reader
        self.runner = _LifecycleRunner(callback)
        return self.runner


class _RecordingManager(ExecutionManager):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.finalizing_seen = False

    def finalize(self, record, **kwargs):
        assert self.store.get(record.id).state == RunState.FINALIZING
        self.finalizing_seen = True
        return super().finalize(record, **kwargs)


def _live_run(tmp_path, reader):
    output = tmp_path / "run.out"
    output.write_text("captured output")
    store = RunStore(tmp_path / "runs")
    manager = _RecordingManager(store)
    executor = _LifecycleExecutor()
    run = manager.create_run(
        executor,
        reader=reader,
        spec=SimpleNamespace(
            directory=str(tmp_path),
            input_file=None,
            output_file=str(output),
            input_files=(),
            output_files=(),
        ),
        task="DOS",
    )
    reader.store = store
    reader.run_id = run.id
    run.submit()
    return store, manager, executor, run


def _available_record(tmp_path):
    output = tmp_path / "run.out"
    output.write_text("captured output")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_file=str(output), output_files=(str(output),)),
    )
    record = manager.transition(record, state=RunState.FINALIZING)
    record = manager.finalize(record, state=RunState.AVAILABLE, return_code=0)
    return store, manager, record


def test_successful_finalization_checks_only_after_available(tmp_path):
    reader = _ResultReader()
    store, manager, executor, run = _live_run(tmp_path, reader)

    executor.runner.complete()

    record = store.get(run.id)
    assert manager.finalizing_seen
    assert reader.calls == 1
    assert record.state == RunState.CHECKED
    assert record.failure_reason is None
    assert record.finished_at is not None
    assert record.result_files.captured_at is not None


def test_retryable_first_check_stays_available_until_explicit_retry(tmp_path):
    reader = _ResultReader(errors=(OutputParseError("truncated output"),))
    store, manager, executor, run = _live_run(tmp_path, reader)

    with pytest.raises(PersistentResultMismatchError, match="truncated output"):
        executor.runner.complete()

    available = store.get(run.id)
    assert manager.finalizing_seen
    assert reader.calls == 1
    assert available.state == RunState.AVAILABLE
    assert "truncated output" in available.failure_reason

    # Nothing retries merely because the record is inspected/reloaded.
    assert store.get(run.id).state == RunState.AVAILABLE
    assert reader.calls == 1

    result = load_persistent_result(manager, available, reader=reader)
    checked = store.get(run.id)
    assert reader.calls == 2
    assert result.run_record.state == RunState.CHECKED
    assert checked.state == RunState.CHECKED
    assert checked.failure_reason is None


def test_unexpected_first_check_failure_is_terminal_failed(tmp_path):
    error = RuntimeError("reader bug")
    reader = _ResultReader(errors=(error,))
    store, _manager, executor, run = _live_run(tmp_path, reader)

    with pytest.raises(RuntimeError, match="reader bug") as exc_info:
        executor.runner.complete()

    assert exc_info.value is error
    failed = store.get(run.id)
    assert failed.state == RunState.FAILED
    assert failed.failure_reason == "reader bug"
    assert reader.calls == 1


def test_failed_execution_is_not_result_checked(tmp_path):
    reader = _ResultReader()
    store, manager, executor, run = _live_run(tmp_path, reader)

    executor.runner.complete(RunState.FAILED)

    record = store.get(run.id)
    assert manager.finalizing_seen
    assert record.state == RunState.FAILED
    assert record.return_code == 1
    assert reader.calls == 0


def test_post_backend_finalization_failure_does_not_check_result(tmp_path):
    reader = _ResultReader()
    store, manager, executor, run = _live_run(tmp_path, reader)

    executor.runner.complete(error=OSError("result directory disappeared"))

    record = store.get(run.id)
    assert manager.finalizing_seen
    assert record.state == RunState.FINALIZATION_FAILED
    assert record.failure_reason == "result directory disappeared"
    assert record.result_files.captured_at is None
    assert reader.calls == 0


def test_successful_finalization_clears_stale_failure_reason(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(
        record,
        state=RunState.FINALIZING,
        message="temporary backend message",
    )

    record = manager.finalize(record, state=RunState.AVAILABLE, return_code=0)

    assert record.state == RunState.AVAILABLE
    assert record.failure_reason is None


def test_retryable_historical_reconstruction_error_stays_available(tmp_path):
    output = tmp_path / "run.out"
    output.write_text("captured output")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="NO_SUCH_TASK",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_file=str(output), output_files=(str(output),)),
    )
    record = manager.transition(record, state=RunState.FINALIZING)
    record = manager.finalize(record, state=RunState.AVAILABLE, return_code=0)

    with pytest.raises(UnrecoverablePersistentResultError, match="unknown persisted task"):
        check_persistent_result(manager, record)

    current = store.get(record.id)
    assert current.state == RunState.AVAILABLE
    assert "unknown persisted task" in current.failure_reason
    assert record.id not in {item.id for item in store.runs(active=True)}


def test_successful_on_demand_check_marks_checked_and_clears_error(tmp_path):
    store, manager, record = _available_record(tmp_path)
    record = manager.record_result_check_failure(
        record,
        PersistentResultUnavailableError(record, "temporary parser problem"),
        retryable=True,
    )
    reader = _ResultReader(store)
    reader.run_id = record.id

    result = check_persistent_result(manager, record, reader=reader)

    current = store.get(record.id)
    assert result.run_record.state == RunState.CHECKED
    assert current.state == RunState.CHECKED
    assert current.failure_reason is None


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


def test_legacy_finished_record_is_upgraded_and_checked_on_demand(tmp_path):
    output = tmp_path / "run.out"
    output.write_text("legacy output")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_file=str(output), output_files=(str(output),)),
    )
    record = manager.transition(record, state=RunState.FINISHED, return_code=0)
    reader = _ResultReader(store)
    reader.run_id = record.id

    result = load_persistent_result(manager, record, reader=reader)

    current = store.get(record.id)
    assert result.run_record.state == RunState.CHECKED
    assert current.state == RunState.CHECKED
    assert reader.calls == 1
