from types import SimpleNamespace

import pytest
from kkrunner import RunState, RunStore
from kkrunner.files import FileManifest
from kkrunner.runs.manager import RunManager

from ase2sprkkr.outputs.project_results import ProjectResultUnavailableError, load_project_result


def test_later_primary_output_is_not_mistaken_for_uncaptured_history(tmp_path):
    store = RunStore(tmp_path / 'runs')
    manager = RunManager(store)
    input_file = tmp_path / 'run.inp'
    output_file = tmp_path / 'run.out'
    input_file.write_text('input')
    spec = SimpleNamespace(directory=str(tmp_path), input_file=str(input_file), output_file=str(output_file), input_files=(str(input_file),), output_files=(str(output_file),), cache_files=())
    record = manager.new_record(directory=tmp_path, files=FileManifest.from_spec(spec), metadata={'ase2sprkkr': {'task': str('DOS').upper()}}, processing_required=True)
    store.save(record)
    record = manager.complete(record, state=RunState.FINISHED, return_code=0, artifacts_available=False)
    output_file.write_text('output from a later calculation')

    class Reader:

        def read_from_file(self, *args, **kwargs):
            raise AssertionError('uncaptured output must never be parsed')
    with pytest.raises(ProjectResultUnavailableError, match='unavailable'):
        load_project_result(manager, record.id)
