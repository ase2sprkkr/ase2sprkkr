import os

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters


@pytest.mark.parametrize("task", ["BSF", "BSFEK", "BSFKK"])
def test_bsf_aliases_share_native_output_manifest(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"

    names = {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }

    assert names == {
        "manifest_test_BLOCHSF.bsf",
        "manifest_test_BLOCHSF_spol.bsf",
    }
