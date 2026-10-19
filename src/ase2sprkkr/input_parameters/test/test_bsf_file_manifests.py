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

def test_bsf_schema_fallback_preserves_native_task_selection():
    parameters = InputParameters.create("BSF")
    assert parameters.definition_name == "BSF"
    assert parameters.task_name == "BSFKK"
    assert parameters.task_definition.name == "BSFKK"

    parameters.ENERGY.NE = [200]
    assert parameters.task_name == "BSFEK"
    assert InputParameters.create_input_parameters("BSF").definition_name == "BSF"


def test_unknown_task_or_definition_still_fails():
    from ase2sprkkr.common.module_registry import RegistryEntryNotFound

    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create("DOES_NOT_EXIST")
