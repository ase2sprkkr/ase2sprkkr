from ase2sprkkr.input_parameters.input_parameters import InputParameters


def test_dos_output_is_namespaced_by_dataset(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe_DOS.dos"),)


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

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe_XAS.tau.dir"),)


def test_bsfek_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("BSFEK")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe_BLOCHSF.tau.dir"),)


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
        str(tmp_path / "Fe_JXC_XCPLTEN_Jij.dat"),
        str(tmp_path / "Fe_JXC_XCPLTEN_Jij_*.agr"),
        str(tmp_path / "Fe_JXC_XCPLTEN_Dij.dat"),
        str(tmp_path / "Fe_JXC_XCPLTEN_Dij_[xyz]_*.agr"),
    )


def test_jxc_sp_srel_declares_scalar_and_plot_outputs(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.MODE.MODE = "SP-SREL"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe_JXC_J_ij.dat"),
        str(tmp_path / "Fe_JXC_Jij_*.agr"),
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
        outputs = (
            str(tmp_path / "Fe_SPEC.out"),
            str(tmp_path / "Fe_data.spc"),
        )
        if task == "CLPES":
            outputs += (str(tmp_path / "Fe_corewf"),)
        assert parameters.output_files(str(tmp_path)) == outputs

        # Selecting S polarization does not itself enable two output channels.
        parameters.SPEC_PH.POL_P = "S"
        assert parameters.output_files(str(tmp_path)) == outputs

        if task != "CLPES":
            for npol in (0, 3):
                parameters.SPEC_PH.NPOL = npol
                output_files = parameters.output_files(str(tmp_path))
                expected = (
                    str(tmp_path / "Fe_SPEC.out"),
                    str(tmp_path / "Fe_data_pol1.spc"),
                    str(tmp_path / "Fe_data_pol2.spc"),
                )
                assert output_files == expected

    parameters = InputParameters.create("BAND")
    parameters.CONTROL.DATASET = "Fe"
    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe_SPEC.out"), str(tmp_path / "Fe_data.spc"),
    )

    parameters = InputParameters.create("SPLEED")
    parameters.CONTROL.DATASET = "Fe"
    for diagram in (0, 1, 2):
        parameters.SPEC_EL.TYP = diagram
        declared = parameters.output_files(str(tmp_path))
        assert tuple(map(str, declared)) == (
            str(tmp_path / "Fe_SPEC.out"), str(tmp_path / "Fe_data.spc"),
        )
        assert all(not f.optional for f in declared)


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
    parameters.ENERGY.ImE = 0.0
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

def test_clpes_core_wavefunction_is_a_declared_output(tmp_path):
    parameters = InputParameters.create("CLPES")
    parameters.CONTROL.DATASET = "Fe"
    assert str(tmp_path / "Fe_corewf") in parameters.output_files(str(tmp_path))


def test_native_dataset_and_explicit_tau_paths(tmp_path):
    from ase2sprkkr.input_parameters.task_definition import format_file_template

    parameters = InputParameters.create("XAS")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = "CUSTOM"
    assert format_file_template(parameters, "{dataset}_{adsi}") == "Fe_CUSTOM"
    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "Fe_CUSTOM.tau.dir"),
    )

    parameters.TAU.TAUDIR = "existing.tau.dir"
    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "existing.tau.dir"),
    )
    parameters.TAU.TAUPATH = "shared"
    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "shared" / "existing.tau.dir"),
    )


def test_native_adsi_input_output_and_optional_tau(tmp_path):
    for name in ("AES", "COMPTON", "POSANI", "VBPES"):
        parameters = InputParameters.create(name)
        parameters.CONTROL.DATASET = "Fe"
        parameters.CONTROL.ADSI = "OTHER"
        assert str(tmp_path / "Fe_OTHER.tau.dir") in parameters.cache_files(str(tmp_path))
        assert str(tmp_path / "Fe_OTHER.tau.dir") not in parameters.input_files(str(tmp_path))

    for name in ("APS", "NRAES"):
        parameters = InputParameters.create(name)
        parameters.CONTROL.DATASET = "Fe"
        parameters.CONTROL.ADSI = "OTHER"
        assert parameters.input_files(str(tmp_path)) == (
            str(tmp_path / "Fe_OTHER.dos"),
        )

    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = "OTHER"
    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe_OTHER.dos"),
    )


def test_native_adsi_exchange_and_original_dataset_template(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.ADSI = "CUSTOM"
    parameters.MODE.MODE = "FREL"
    assert str(tmp_path / "Fe_CUSTOM_XCPLTEN_Jij.dat") in parameters.output_files(str(tmp_path))



def test_producer_consumer_native_dataset_defaults(tmp_path):
    producer = InputParameters.create("DOS")
    producer.CONTROL.DATASET = "Fe"
    produced = producer.output_files(str(tmp_path))

    for task in ("APS", "NRAES"):
        consumer = InputParameters.create(task)
        consumer.CONTROL.DATASET = "Fe"
        assert consumer.TASK.TASK() == task
        assert consumer.input_files(str(tmp_path)) == produced
        assert consumer.CONTROL.ADSI() == "DOS"

    producer = InputParameters.create("POSANIPREP")
    consumer = InputParameters.create("POSANI")
    producer.CONTROL.DATASET = "Fe"
    consumer.CONTROL.DATASET = "Fe"
    assert producer.TASK.TASK() == "POSANIPREP"
    assert consumer.TASK.TASK() == "POSANI"
    assert producer.CONTROL.ADSI() == consumer.CONTROL.ADSI() == "POSANI"
    assert producer.output_files(str(tmp_path)) == consumer.input_files(str(tmp_path))


def test_native_tau_absolute_path_and_native_prefix_semantics(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.CONTROL.DATASET = "Fe"
    absolute = str(tmp_path / "saved.tau.dir")
    parameters.TAU.TAUDIR = absolute
    assert parameters.cache_files(str(tmp_path)) == (absolute,)

    parameters.TAU.TAUPATH = "remote"
    # Native concatenates TAUPATH and TAUDIR; it does not discard TAUPATH
    # when TAUDIR happens to start with a slash.
    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "remote" / absolute.lstrip("/")),
    )


def test_native_datset0_remains_available_to_legacy_outputs(tmp_path):
    for task, suffix in (("JALF", "_JALF.dat"), ("JORB", "_JORB.dat")):
        params = InputParameters.create(task)
        params.CONTROL.DATASET = "Fe"
        params.CONTROL.ADSI = "CUSTOM"
        assert str(tmp_path / ("Fe" + suffix)) in params.output_files(str(tmp_path))


def test_native_plotprs_uses_datset0_even_with_adsi(tmp_path):
    # fpplotprs.f calls XMGRHEAD(TRIM(DATSET0), ...).
    for task in ("PLOTPRS", "PLOTPOT", "PLOTRHO", "PLOTSFN"):
        params = InputParameters.create(task)
        params.CONTROL.DATASET = "Fe"
        params.CONTROL.ADSI = "CUSTOM"
        for output in params.output_files(str(tmp_path)):
            assert str(output).startswith(str(tmp_path / "Fe_"))
            assert not str(output).startswith(str(tmp_path / "Fe_CUSTOM_"))


def test_bsf_output_declarations_follow_custom_adsi(tmp_path):
    for task in ("BSFEK", "BSFKK"):
        params = InputParameters.create(task)
        params.CONTROL.DATASET = "Fe"
        params.CONTROL.ADSI = "CUSTOM"
        declared = params.output_files(str(tmp_path))
        assert str(declared[0]) == str(tmp_path / "Fe_CUSTOM.bsf")
        assert str(declared[1]) == str(tmp_path / "Fe_CUSTOM_spol.bsf")
        assert declared[1].optional


def test_separate_dataset_adsi_placeholders_and_empty_suffix():
    from types import SimpleNamespace
    from ase2sprkkr.input_parameters.task_definition import format_file_template

    control = SimpleNamespace(DATASET=lambda: "Fe", ADSI=lambda: "CUSTOM")
    params = SimpleNamespace(CONTROL=control)
    assert format_file_template(params, "{dataset}_{adsi}.dos") == "Fe_CUSTOM.dos"
    assert format_file_template(params, "{dataset}.dos") == "Fe.dos"
    assert format_file_template(params, "{adsi}") == "CUSTOM"
    control.ADSI = lambda: ""
    assert format_file_template(params, "{dataset}_{adsi}.dos") == "Fe.dos"
    assert format_file_template(params, "{dataset}") == "Fe"
    assert format_file_template(params, "{adsi}") == ""


def test_file_callbacks_return_literal_paths(tmp_path):
    from ase2sprkkr.input_parameters.task_definition import DeclaredFile

    params = InputParameters.create("DOS")
    params.CONTROL.DATASET = "Fe"
    params.CONTROL.ADSI = "DOS"
    resolve = InputParameters.input_parameters_registry._resolve_file_templates
    results = resolve(params, str(tmp_path), (
        "{dataset}_{adsi}.dos",
        lambda _: DeclaredFile("archive_{2026}.dat", optional=True),
        lambda _: "report_{custom}.dat",
    ))
    assert results == (
        str(tmp_path / "Fe_DOS.dos"),
        DeclaredFile(str(tmp_path / "archive_{2026}.dat"), optional=True),
        str(tmp_path / "report_{custom}.dat"),
    )
    assert results[1].optional
