import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(parameters, tmp_path, kind="output"):
    return {
        os.path.basename(os.fspath(path))
        for path in getattr(parameters, f"{kind}_files")(tmp_path)
    }


def test_gilbert_tau_cache_and_output_families(tmp_path):
    parameters = InputParameters.create("GILBERT")
    parameters.CONTROL.DATASET = "manifest_test"
    assert "manifest_test.tau.dir" in _names(parameters, tmp_path, "cache")
    assert {
        "manifest_testalfa_x.dat",
        "manifest_testalfa_T.dat",
        "manifest_test_alfa_T.agr",
    } <= _names(parameters, tmp_path)


def test_stoner_diagnostic_output_families(tmp_path):
    parameters = InputParameters.create("STONER")
    parameters.CONTROL.DATASET = "manifest_test"
    assert {
        "manifest_test_*_g-spin.dat",
        "manifest_test_*_g-orb.dat",
        "manifest_test_*_VWN_A-xc.dat",
    } <= _names(parameters, tmp_path)
