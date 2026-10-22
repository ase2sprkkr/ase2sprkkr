import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"
    return {
        os.path.basename(os.fspath(path))
        for path in parameters.output_files(tmp_path)
    }


def test_fscat_output_family(tmp_path):
    assert "manifest_test_GEO?_*.agr" in _names("FSCAT", tmp_path)


def test_nucsscpl_outputs(tmp_path):
    assert {
        "manifest_test_A_ij.dat",
        "RK.dat",
        "tab_aij.dat",
        "manifest_test_Aij_*.agr",
        "manifest_test_Rij_*.agr",
    } <= _names("NUCSSCPL", tmp_path)
