import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"
    return {
        os.path.basename(os.fspath(path))
        for path in parameters.output_files(tmp_path)
    }


def test_optabs_output_families(tmp_path):
    assert {
        "sig_?????_ISPR??_NV.dat",
        "sig_?????_ISPR??_VC.dat",
    } <= _names("OPTABS", tmp_path)


def test_bsdosij_output_families(tmp_path):
    assert {
        "manifest_test_BSDOSIJ_XCPLTEN_Jij.dat",
        "manifest_test_BSDOSIJ_XCPLTEN_Jij_*.agr",
        "manifest_test_BSDOSIJ_XCPLTEN_Dij_[xyz]_*.agr",
    } <= _names("BSDOSIJ", tmp_path)
