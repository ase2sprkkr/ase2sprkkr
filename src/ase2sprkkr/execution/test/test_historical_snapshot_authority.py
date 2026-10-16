from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager,
    FileState,
    ResultFilesRecord,
    RunState,
    RunStore,
)
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    load_persistent_result,
)


def test_later_primary_output_is_not_mistaken_for_uncaptured_history(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("input")
    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
        cache_files=(),
    )
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FINISHED, return_code=0)

    # The original run produced no output. Reusing the directory later must not
    # make this unrelated file part of the historical result.
    output_file.write_text("output from a later calculation")

    class Reader:
        def read_from_file(self, *args, **kwargs):
            raise AssertionError("uncaptured output must never be parsed")

    with pytest.raises(PersistentResultUnavailableError, match="unavailable"):
        load_persistent_result(manager, record.id, reader=Reader())


def test_later_cache_update_does_not_overwrite_historical_result(tmp_path):
    cache = tmp_path / "shared.tau"
    output = tmp_path / "result.dat"
    cache.write_text("cache from this run")
    output.write_text("immutable result")

    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    files = ResultFilesRecord(
        output_files=(str(cache), str(output)),
        cache_files=(str(cache),),
    )
    record = manager.new_record(
        task="BSFEK",
        directory=tmp_path,
        result_files=files,
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FINISHED, return_code=0)

    assert {version.path for version in record.result_files.output_versions} == {
        str(output)
    }

    cache.write_text("cache updated by a later run")
    assert record.result_files.state is FileState.CURRENT
