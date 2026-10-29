import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def test_chi_standard_outputs(tmp_path):
    parameters = InputParameters.create("CHI")
    parameters.CONTROL.DATASET = "manifest_test"
    parameters.CONTROL.ADSI = ""
    names = {
        os.path.basename(os.fspath(path))
        for path in parameters.output_files(tmp_path)
    }
    assert {
        "manifest_test*_dens.agr",
        "manifest_test*_chifmag.agr",
        "manifest_test_*_GAMMA-OP.dat",
        "manifest_test_*_CHI-r.dat",
        "manifest_test_*_*_I-spin.dat",
    } <= names
