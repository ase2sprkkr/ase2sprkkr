from ase2sprkkr.input_parameters.input_parameters import InputParameters


def test_dos_output_is_namespaced_by_dataset(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.dos"),)


def test_scf_declares_new_potential_as_output(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = "Fe.pot"

    assert tuple(map(str, parameters.output_files(str(tmp_path)))) == (
        str(tmp_path / "Fe.pot_new"), str(tmp_path / "*_SCFSTART.dos"),
    )


def test_scf_fullpot_shape_function_is_optional_cache(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.SCF.FULLPOT = True
    parameters.CONTROL.SFNFIL = "Fe.sfn"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.sfn"),)


def test_xas_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.tau.dir"),)


def test_bsfek_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("BSFEK")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.tau.dir"),)


def test_gilbert_fluctuation_file_is_required_for_m_t(tmp_path):
    parameters = InputParameters.create("GILBERT")
    parameters.TASK.SETFLUCT = "M_T"
    parameters.TASK.FLUCTFIL = "magnetization.dat"

    assert parameters.input_files(str(tmp_path)) == (
        str(tmp_path / "magnetization.dat"),
    )


def test_phagen_outputs_follow_explicit_potential_name(tmp_path):
    parameters = InputParameters.create("PHAGEN")
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.CONTROL.DATASET = "unrelated"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe.pot_ex_sprkkr"),
        str(tmp_path / "Fe.pot.only_coulomb"),
    )


def test_jxc_frel_declares_tensor_and_plot_outputs(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.MODE.MODE = "FREL"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe_XCPLTEN_Jij.dat"),
        str(tmp_path / "Fe_XCPLTEN_Jij_*.agr"),
        str(tmp_path / "Fe_XCPLTEN_Dij.dat"),
        str(tmp_path / "Fe_XCPLTEN_Dij_[xyz]_*.agr"),
    )


def test_jxc_sp_srel_declares_scalar_and_plot_outputs(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.MODE.MODE = "SP-SREL"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe_J_ij.dat"),
        str(tmp_path / "Fe_Jij_*.agr"),
    )


def test_spec_inputs_follow_native_flags(tmp_path):
    for task in ("ARPES", "AIPES", "BAND", "SPLEED", "CLPES"):
        parameters = InputParameters.create(task)
        parameters.CONTROL.DATASET = "Fe"
        assert parameters.input_files(str(tmp_path)) == ()
        parameters.TASK.INPVER = 0
        assert parameters.input_files(str(tmp_path)) == (
            str(tmp_path / "input.inp"),
        )
        parameters.TASK.STRVER = 0
        assert parameters.input_files(str(tmp_path)) == (
            str(tmp_path / "input.inp"), str(tmp_path / "struc.inp"),
        )


def test_spec_outputs_follow_polarization_and_task(tmp_path):
    for task in ("ARPES", "AIPES", "CLPES"):
        parameters = InputParameters.create(task)
        parameters.CONTROL.DATASET = "Fe"
        assert parameters.output_files(str(tmp_path)) == (
            str(tmp_path / "Fe_SPEC.out"), str(tmp_path / "Fe_data.spc"),
        )
        parameters.SPEC_PH.POL_P = "S"
        assert parameters.output_files(str(tmp_path)) == (
            str(tmp_path / "Fe_SPEC.out"),
            str(tmp_path / "Fe_data_pol1.spc"),
            str(tmp_path / "Fe_data_pol2.spc"),
        )

    for task in ("BAND", "SPLEED"):
        parameters = InputParameters.create(task)
        parameters.CONTROL.DATASET = "Fe"
        assert parameters.output_files(str(tmp_path)) == (
            str(tmp_path / "Fe_SPEC.out"), str(tmp_path / "Fe_data.spc"),
        )

def test_scf_invalid_grid_and_mixing_are_rejected():
    import pytest
    from ase2sprkkr.common.warnings import DataValidityError

    parameters = InputParameters.create("SCF")
    parameters.ENERGY.GRID = 3
    parameters.ENERGY.NE = 1
    with pytest.raises(DataValidityError, match="NE>=2"):
        parameters.validate("save")

    parameters = InputParameters.create("SCF")
    parameters.SCF.FULLPOT = True
    parameters.SCF.ALG = "TCHEBY"
    with pytest.raises(DataValidityError, match="TCHEBY"):
        parameters.validate("save")

def test_dos_native_energy_constraints():
    import pytest
    from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning

    parameters = InputParameters.create("DOS")
    parameters.CONTROL.POTFIL = "Fe.pot"
    parameters.ENERGY.NE = 1
    with pytest.raises(DataValidityError, match="two energy points"):
        parameters.validate("save")
    parameters.ENERGY.NE = 300
    parameters.ENERGY.ImE = 0
    with pytest.warns(DataValidityWarning, match="real axis"):
        parameters.validate("save")


def test_scf_startup_dos_is_optional(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = "Fe.pot"
    files = parameters.output_files(str(tmp_path))
    assert next(f for f in files if str(f).endswith("_SCFSTART.dos")).optional
    assert not next(f for f in files if str(f).endswith("Fe.pot_new")).optional


def test_bsf_output_file_declarations(tmp_path):
    for task in ("BSFEK", "BSFKK"):
        parameters = InputParameters.create(task)
        parameters.CONTROL.DATASET = "Fe"
        assert tuple(map(str, parameters.output_files(str(tmp_path)))) == (
            str(tmp_path / "Fe_BLOCHSF.bsf"),
            str(tmp_path / "Fe_BLOCHSF_spol.bsf"),
        )
