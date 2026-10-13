from ase2sprkkr.input_parameters import InputParameters


def test_dos_output_is_namespaced_by_dataset(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.dos"),)


def test_scf_declares_new_potential_as_output(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = "Fe.pot"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.pot_new"),)


def test_scf_fullpot_shape_function_is_optional_cache(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.SCF.FULLPOT = True
    parameters.CONTROL.SFNFIL = "Fe.sfn"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.sfn"),)


def test_xas_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.tau.dir"),)
