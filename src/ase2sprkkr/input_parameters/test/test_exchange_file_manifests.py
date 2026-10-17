import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


DATASET = "manifest_test"


def _names(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = DATASET
    return parameters, {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }


def test_jalf_and_jorb_native_names(tmp_path):
    _, jalf = _names("JALF", tmp_path)
    assert {f"{DATASET}_JALF.dat", "curdns1", "curdns2", "curdns3"} <= jalf
    assert f"{DATASET}0_JALF.dat" not in jalf

    _, jorb = _names("JORB", tmp_path)
    assert {
        f"{DATASET}_JORB.dat",
        f"{DATASET}_JORB_UP.dat",
        f"{DATASET}_JORB_DN.dat",
    } <= jorb
    assert f"{DATASET}0_JORB.dat" not in jorb


def test_jxc_frel_outputs(tmp_path):
    _, names = _names("JXC", tmp_path)
    assert {
        f"{DATASET}_XCPLTEN_Jij.dat",
        f"{DATASET}_XCPLTEN_Jij_*.agr",
        f"{DATASET}_XCPLTEN_Dij.dat",
        f"{DATASET}_XCPLTEN_Dij_[xyz]_*.agr",
    } <= names


def test_jxc_sp_srel_outputs(tmp_path):
    parameters, _ = _names("JXC", tmp_path)
    parameters.MODE.MODE = "SP-SREL"
    names = {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }
    assert {f"{DATASET}_J_ij.dat", f"{DATASET}_Jij_*.agr"} <= names
    assert f"{DATASET}_XCPLTEN_Jij.dat" not in names


def test_dmi_outputs(tmp_path):
    _, names = _names("DMI", tmp_path)
    assert {
        f"{DATASET}_DMIVEC_Dij.dat",
        f"{DATASET}_DMIVEC_Dij_[xyz]_*.agr",
    } <= names
