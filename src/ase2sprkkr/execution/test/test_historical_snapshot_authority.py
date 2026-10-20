from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import ExecutionSpec, RunState, RunStore, FileManifest
from ase2sprkkr.execution.manager import ExecutionManager
from ase2sprkkr.outputs.project_results import (
    ProjectResultUnavailableError,
    load_project_result,
)


def test_later_primary_output_is_not_mistaken_for_uncaptured_history(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("input")
    spec = ExecutionSpec(
        argv=(),
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
    )
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        files=FileManifest.from_spec(spec),
    )
    store.save(record)
    record = manager.complete(record, state=RunState.FINISHED, return_code=0, artifacts_available=False)

    # The original run produced no output. Reusing the directory later must not
    # make this unrelated file part of the historical result.
    output_file.write_text("output from a later calculation")

    class Reader:
        def read_from_file(self, *args, **kwargs):
            raise AssertionError("uncaptured output must never be parsed")

    with pytest.raises(ProjectResultUnavailableError, match="unavailable"):
        load_project_result(manager, record.id)
