import pickle

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.input_parameters_definitions import InputParametersDefinition
from ase2sprkkr.input_parameters.task_definition import TaskFile


def test_task_file_rejects_empty_type_identifier():
    with pytest.raises(ValueError, match="type must not be empty"):
        TaskFile("result.dat", type="   ")


def test_input_parameter_schemas_do_not_own_task_file_manifests():
    for name in ("input_files", "output_files", "cache_files"):
        assert not hasattr(InputParametersDefinition, name)


def test_builtin_task_file_metadata_is_pickleable():
    for task in InputParameters.task_registry.values():
        restored = pickle.loads(pickle.dumps(task))
        assert restored == task
