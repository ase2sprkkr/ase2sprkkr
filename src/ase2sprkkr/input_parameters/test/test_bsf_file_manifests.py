import os

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters


@pytest.mark.parametrize("task", ["BSFEK", "BSFKK"])
def test_bsf_aliases_share_native_output_manifest(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"

    names = {
        os.path.basename(os.fspath(path))
        for path in parameters.output_files(tmp_path)
    }

    assert names == {
        "manifest_test_BLOCHSF.bsf",
        "manifest_test_BLOCHSF_spol.bsf",
    }

def test_bsf_schema_access_is_distinct_from_task_creation():
    from ase2sprkkr.common.module_registry import RegistryEntryNotFound

    # The bare BSF name identifies an input schema, not an executable task.
    schema = InputParameters.definition("BSF")
    parameters = schema.create_object()
    assert parameters.definition_name == "BSF"
    assert parameters.task_name == "BSFKK"
    assert parameters.task_definition.name == "BSFKK"
    parameters.ENERGY.NE = [200]
    assert parameters.task_name == "BSFEK"
    assert parameters.task_definition.name == "BSFEK"

    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create("BSF")
    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create_input_parameters("BSF")


def test_unknown_task_or_definition_still_fails():
    from ase2sprkkr.common.module_registry import RegistryEntryNotFound

    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create("DOES_NOT_EXIST")
