from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ArtifactAvailabilityState,
    ExecutionManager,
    ExecutionSpec,
    Executor,
    ExecutorPool,
    FileState,
    RunReference,
    RunState,
    RunStore,
)
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.outputs import persistent as persistent_results
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    load_persistent_result,
)


class _HistoryExecutor(Executor):
    name = "history"

    def __init__(self):
        self.fetch_calls = []
        self.output_file = None
        self.restored_output = "original output"

    def fetch_result_files(self, reference, files):
        files = tuple(files)
        self.fetch_calls.append((reference, files))
        if self.output_file is None or self.output_file not in files:
            return
        with open(self.output_file, "w") as output:
            output.write(self.restored_output)


def _spec(tmp_path):
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("input")
    return ExecutionSpec(
        argv=(),
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
    )


def _manager(store):
    return ExecutionManager(store, executor_pool=ExecutorPool())


def _terminal_run(tmp_path, *, state=RunState.FINISHED, return_code=0, task="DOS"):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    executor = _HistoryExecutor()
    spec = _spec(tmp_path)
    output_file = tmp_path / "run.out"
    output_file.write_text("original output")
    executor.output_file = str(output_file)

    endpoint = endpoint_for_executor(executor)
    manager.register_executor(executor, endpoint)
    reference = RunReference(
        executor="history",
        native_id="job-1",
        directory=str(tmp_path),
        output_file=str(output_file),
    )
    record = manager.new_record(
        task=task,
        directory=tmp_path,
        endpoint=endpoint,
        result_files=manager._result_files_for_spec(spec),
    )
    record = manager.store.save(record)
    record = manager.complete(
        record,
        state=state,
        reference=reference,
        return_code=return_code,
        artifacts_available=True,
    )
    return manager, executor, record


def _patch_reader(monkeypatch, parsed):
    class Reader:
        input_parameters = object()

        def read_from_file(self, output, **kwargs):
            parsed.append((output, kwargs))
            return SimpleNamespace(return_code=kwargs["return_code"])

    monkeypatch.setattr(
        persistent_results,
        "_reader_for_record",
        lambda manager, record: Reader(),
    )


def test_history_listing_does_not_fetch_missing_result_files(tmp_path):
    manager, executor, record = _terminal_run(tmp_path)
    (tmp_path / "run.out").unlink()

    records = manager.store.runs()

    assert records[0].id == record.id
    assert executor.fetch_calls == []


def test_load_persistent_result_restores_exact_missing_historical_output(
    tmp_path, monkeypatch
):
    manager, executor, record = _terminal_run(tmp_path)
    (tmp_path / "run.out").unlink()
    parsed = []
    _patch_reader(monkeypatch, parsed)

    result = load_persistent_result(manager, record.id)

    assert executor.fetch_calls == [
        (record.reference, (str(tmp_path / "run.out"),))
    ]
    assert (tmp_path / "run.out").read_text() == "original output"
    assert parsed[0][0] == str(tmp_path / "run.out")
    assert result.return_code == 0
    assert result.run_record.state == RunState.FINISHED
    assert result.run_record.result_files.state == FileState.CURRENT
    assert (
        result.run_record.result_files.availability_state
        == ArtifactAvailabilityState.AVAILABLE
    )


def test_restored_changed_content_is_rejected_as_overwritten(tmp_path, monkeypatch):
    manager, executor, record = _terminal_run(tmp_path)
    (tmp_path / "run.out").unlink()
    executor.restored_output = "changed output with different contents"

    class Reader:
        input_parameters = object()

        def read_from_file(self, *args, **kwargs):
            raise AssertionError("overwritten historical output must not be parsed")

    monkeypatch.setattr(
        persistent_results,
        "_reader_for_record",
        lambda manager, record: Reader(),
    )

    with pytest.raises(PersistentResultUnavailableError, match="overwritten"):
        load_persistent_result(manager, record.id)

    stored = manager.store.get(record.id)
    assert stored.state == RunState.FINISHED
    assert stored.result_files.state == FileState.OVERWRITTEN


def test_retrieval_programming_error_is_not_persisted_as_artifact_failure(tmp_path):
    manager, executor, record = _terminal_run(tmp_path)
    (tmp_path / "run.out").unlink()

    def broken_fetch(reference, files):
        raise KeyError("programming error")

    executor.fetch_result_files = broken_fetch

    with pytest.raises(KeyError, match="programming error"):
        manager.retrieve_result_files(record.id)

    stored = manager.store.get(record.id)
    assert (
        stored.result_files.availability_state
        is ArtifactAvailabilityState.AVAILABLE
    )


def test_finished_history_without_available_output_is_not_synthesized(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = _manager(store)
    spec = _spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    record = manager.store.save(record)
    record = manager.complete(
        record,
        state=RunState.FINISHED,
        return_code=0,
        artifacts_available=False,
    )

    with pytest.raises(PersistentResultUnavailableError, match="unavailable"):
        load_persistent_result(manager, record.id)


def test_failed_history_without_return_code_is_not_parsed_as_success(
    tmp_path, monkeypatch
):
    manager, _, record = _terminal_run(
        tmp_path,
        state=RunState.FAILED,
        return_code=None,
    )
    parsed = []
    _patch_reader(monkeypatch, parsed)

    result = load_persistent_result(manager, record.id)

    assert parsed[0][1]["return_code"] != 0
    assert result.return_code != 0
    assert result.run_record.state == RunState.FAILED
