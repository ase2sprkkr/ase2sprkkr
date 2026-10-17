import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(parameters, tmp_path):
    return {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }


def test_chixas_output_families(tmp_path):
    parameters = InputParameters.create("CHIXAS")
    parameters.CONTROL.DATASET = "manifest_test"
    names = _names(parameters, tmp_path)
    assert "manifest_test_*_L23.rat" in names
    assert {f"manifest_test_*_{kind}.pol" for kind in ("ss", "ls", "sl", "ll")} <= names
    assert {f"manifest_test_*_L23_spec_{kind}" for kind in ("0", "s", "o", "sr")} <= names
    assert {
        f"manifest_test_*_L23_{kind}.rat"
        for kind in ("SSZ", "SSX", "SST", "LLZ", "LLX", "LLT", "TOT")
    } <= names
    assert "manifest_test_*_GG.dat" not in names


def test_chixas_shell_and_print_dependent_outputs(tmp_path):
    parameters = InputParameters.create("CHIXAS")
    parameters.CONTROL.DATASET = "manifest_test"
    parameters.TASK.CL = "3D"
    parameters.TASK.PRINT = 1
    names = _names(parameters, tmp_path)
    assert "manifest_test_*_M45.rat" in names
    assert "manifest_test_*_GG.dat" in names
