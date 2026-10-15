import os
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    ExecutorPool,
    FileRetrievalState,
    ResultFilesState,
    RunReference,
    RunState,
    RunStore,
)
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    load_persistent_result,
)
from ase2sprkkr.outputs.task_result import TaskResult


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
        self.exception = None

    def submit(self):
        self.reference = RunReference(executor="history", native_id="job-1")
        self.state = RunState.QUEUED

    def wait(self):
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def finish(self):
        self.result.return_code = 0
        self.state = RunState.FINISHED
        self.callback(self.result)


class _HistoryExecutor(Executor):
    name = "history"

    def __init__(self):
        self.runner = None
        self.fetch_calls = []
        self.output_file = None
        self.restored_output = "original output"
        self.restored_mtime_ns = None

    def create_run(self, reader, spec, callback=None):
        self.output_file = spec.output_file
        self.runner = _Runner(callback)
        return self.runner

    def fetch_result_files(self, reference, files):
        files = tuple(files)
        self.fetch_calls.append((reference, files))
        if self.output_file is None or self.output_file not in files:
            return
        with open(self.output_file, "w") as output:
            output.write(self.restored_output)
        if self.restored_mtime_ns is not None:
            stat = os.stat(self.output_file)
            os.utime(
                self.output_file,
                ns=(stat.st_atime_ns, self.restored_mtime_ns),
            )


class _NoRetrievalExecutor(Executor):
    name = "no-retrieval"


def _spec(tmp_path):
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("input")
    return SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
        cache_files=(),
    )


def _manager(store):
    return ExecutionManager(store, executor_pool=ExecutorPool())


def _finished_run(tmp_path, *, task="DOS"):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    executor = _HistoryExecutor()
    spec = _spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task=task)
    run.submit()
    output_file = tmp_path / "run.out"
    output_file.write_text("original output")
    executor.restored_mtime_ns = output_file.stat().st_mtime_ns
    executor.runner.finish()
    return manager, executor, run


def test_history_listing_does_not_fetch_missing_result_files(tmp_path):
    manager, executor, run = _finished_run(tmp_path)
    (tmp_path / "run.out").unlink()

    records = manager.store.runs()

    assert records[0].id == run.id
    assert executor.fetch_calls == []


def test_load_persistent_result_restores_exact_missing_historical_output(tmp_path):
    manager, executor, run = _finished_run(tmp_path)
    (tmp_path / "run.out").unlink()

    parsed = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            parsed.append((output, kwargs))
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = load_persistent_result(manager, run.id, reader=Reader())

    assert executor.fetch_calls == [
        (run.reference, (str(tmp_path / "run.out"),))
    ]
    assert (tmp_path / "run.out").read_text() == "original output"
    assert parsed[0][0] == str(tmp_path / "run.out")
    assert result.return_code == 0
    assert result.run_record.state == RunState.FINISHED
    assert result.run_record.result_files.state == ResultFilesState.CURRENT
    assert result.run_record.result_files.retrieval_state == FileRetrievalState.AVAILABLE


def test_loaded_persistent_result_can_start_continuation_without_calculator(tmp_path):
    manager, _, run = _finished_run(tmp_path, task="SCF")
    input_parameters = SimpleNamespace(
        CONTROL=SimpleNamespace(POTFIL=lambda: "Fe.pot")
    )

    class Reader:
        def read_from_file(self, output, **kwargs):
            result = TaskResult(
                input_parameters,
                None,
                tmp_path,
                output_file=output,
                task_name="SCF",
            )
            result.complete(None, kwargs["return_code"])
            return result

    result = load_persistent_result(manager, run.id, reader=Reader())

    assert result._calculator is None
    continuation = result.new_task("DOS")
    assert continuation._potential == os.path.join(os.fspath(tmp_path), "Fe.pot")
    assert continuation.input_parameters.task_definition.name == "DOS"


def test_restored_changed_content_is_rejected_as_overwritten(tmp_path):
    manager, executor, run = _finished_run(tmp_path)
    (tmp_path / "run.out").unlink()
    executor.restored_output = "changed output with different size"

    class Reader:
        def read_from_file(self, *args, **kwargs):
            raise AssertionError("overwritten historical output must not be parsed")

    with pytest.raises(PersistentResultUnavailableError, match="overwritten"):
        load_persistent_result(manager, run.id, reader=Reader())

    stored = manager.store.get(run.id)
    assert stored.state == RunState.FINISHED
    assert stored.result_files.state == ResultFilesState.OVERWRITTEN


def test_retrieval_programming_error_is_not_persisted_as_backend_failure(tmp_path):
    manager, executor, run = _finished_run(tmp_path)
    (tmp_path / "run.out").unlink()

    def broken_fetch(reference, files):
        raise KeyError("programming error")

    executor.fetch_result_files = broken_fetch

    with pytest.raises(KeyError, match="programming error"):
        manager.retrieve_result_files(run.id)

    stored = manager.store.get(run.id)
    assert stored.result_files.retrieval_state is not FileRetrievalState.FAILED


def test_finished_history_without_available_output_is_not_synthesized(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    spec = _spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FINISHED, return_code=0)

    class Reader:
        def create_result(self, *args):
            raise AssertionError("successful missing output must not create a shell result")

    with pytest.raises(PersistentResultUnavailableError, match="unavailable"):
        load_persistent_result(manager, record.id, reader=Reader())


def test_failed_history_without_return_code_is_not_parsed_as_success(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    spec = _spec(tmp_path)
    (tmp_path / "run.out").write_text("failed output")
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FAILED)
    parsed = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            parsed.append(kwargs)
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = load_persistent_result(manager, record.id, reader=Reader())

    assert parsed[0]["return_code"] != 0
    assert result.return_code != 0
    assert result.run_record.state == RunState.FAILED


def test_missing_failed_history_uses_inferred_return_code(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    spec = _spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FAILED)

    class Result:
        return_code = None

        def complete(self, error, return_code):
            self.return_code = return_code

    class Reader:
        def create_result(self, output, input_file):
            return Result()

    result = load_persistent_result(manager, record.id, reader=Reader())

    assert result.return_code == 1
    assert result.run_record.state == RunState.FAILED
