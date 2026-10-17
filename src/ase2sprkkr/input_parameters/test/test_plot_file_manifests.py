import os

from ase2sprkkr.input_parameters.input_parameters import InputParameters


def _names(task, tmp_path):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "manifest_test"
    return {
        os.path.basename(os.fspath(path))
        for path in parameters._definition.output_files(parameters, tmp_path)
    }


def test_plotprs_output_families(tmp_path):
    assert {
        "manifest_test_V_*.agr",
        "manifest_test_VNS_*.agr",
        "manifest_test_RHO_*.agr",
        "manifest_test_RNS_*.agr",
        "manifest_test_SFN_M*.agr",
    } <= _names("PLOTPRS", tmp_path)


def test_plot2dpr_output_families(tmp_path):
    assert {
        "PLY-V-ASA-*.ply",
        "POT-ASA-LINES.ply",
        "PLY-V-FP-*.ply",
        "PLY-V-MT-*.ply",
        "PLY-FP-RWS-*.ply",
        "PLY-FP-RATM-*.ply",
        "PLY-V-IS.ply",
        "PLY-FP-CELL-*.ply",
        "PLY-FP-RMT.ply",
        "PLY-FP-CELL-BOUND-V.ply",
        "PLY-FP-CELL-BOUND-0.ply",
        "PLY-FP-RATM.ply",
        "PLY-FP-RWS.ply",
    } <= _names("PLOT2DPR", tmp_path)
