import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"
    return {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }


def test_radial_plot_selectors_have_native_outputs(tmp_path):
    assert {
        "manifest_test_V_*.agr",
        "manifest_test_VNS_*.agr",
    } <= _names("PLOTPOT", tmp_path)
    assert {
        "manifest_test_RHO_*.agr",
        "manifest_test_RNS_*.agr",
    } <= _names("PLOTRHO", tmp_path)
    assert {"manifest_test_SFN_M*.agr"} <= _names("PLOTSFN", tmp_path)


def test_2d_plot_selectors_keep_ply_outputs(tmp_path):
    expected = {
        "PLY-V-ASA-*.ply",
        "POT-ASA-LINES.ply",
        "PLY-V-FP-*.ply",
        "PLY-FP-CELL-*.ply",
    }
    assert expected <= _names("PLOT2DPOT", tmp_path)
    assert expected <= _names("PLOT2DRHO", tmp_path)
