from ..input_parameters import InputParameters


def test_arpes_family_uses_single_spc_output_for_npol_one(tmp_path):
    parameters = InputParameters.create("ARPES")
    parameters.CONTROL.DATASET = "case"
    parameters.SPEC_PH.NPOL = 1

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "case_data.spc"),
    )


def test_arpes_family_uses_polarized_spc_outputs_by_default(tmp_path):
    parameters = InputParameters.create("AIPES")
    parameters.CONTROL.DATASET = "case"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "case_data_pol1.spc"),
        str(tmp_path / "case_data_pol2.spc"),
    )


def test_band_always_uses_single_spc_output(tmp_path):
    parameters = InputParameters.create("BAND")
    parameters.CONTROL.DATASET = "case"
    parameters.SPEC_PH.NPOL = 2

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "case_data.spc"),
    )


def test_arpes_legacy_input_files_follow_inpver_and_strver(tmp_path):
    parameters = InputParameters.create("SPLEED")

    assert parameters.input_files(str(tmp_path)) == ()

    parameters.TASK.INPVER = 0
    assert parameters.input_files(str(tmp_path)) == (
        str(tmp_path / "input.inp"),
        str(tmp_path / "in_bar.inp"),
    )

    parameters.TASK.STRVER = 0
    assert parameters.input_files(str(tmp_path)) == (
        str(tmp_path / "input.inp"),
        str(tmp_path / "in_bar.inp"),
        str(tmp_path / "in_structur.inp"),
    )


def test_torque_declares_tau_as_cache(tmp_path):
    parameters = InputParameters.create("TORQUE")
    parameters.CONTROL.DATASET = "case"

    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "case.tau.dir"),
    )


def test_jxc_declares_optional_dmi_output(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "case"

    assert str(tmp_path / "case_DMIVEC_Dij.dat") not in parameters.output_files(
        str(tmp_path)
    )

    parameters.TASK.DMI = True

    assert str(tmp_path / "case_DMIVEC_Dij.dat") in parameters.output_files(
        str(tmp_path)
    )
