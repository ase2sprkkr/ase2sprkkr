from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.readers.scf import ScfOutputReader
from ase2sprkkr.outputs.task_result import TaskResult


def test_task_potential_production_is_backend_metadata(tmp_path):
    scf = InputParameters.create("SCF")
    dos = InputParameters.create("DOS")

    assert InputParameters.definition("SCF").produces_potential
    assert not InputParameters.definition("DOS").produces_potential
    assert TaskResult(scf, None, tmp_path).produces_potential
    assert not TaskResult(dos, None, tmp_path).produces_potential

    parsed_scf = TaskResult(None, None, tmp_path, task_name="SCF")
    parsed_dos = TaskResult(None, None, tmp_path, task_name="DOS")
    assert parsed_scf.task_name == "SCF"
    assert parsed_scf.produces_potential
    assert not parsed_dos.produces_potential


def test_reader_preserves_task_identity_without_input_parameters(tmp_path):
    reader = ScfOutputReader(None, None, tmp_path, task_name="SCF")
    result = reader.create_result(None)

    assert result.task_name == "SCF"
    assert result.produces_potential


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


def test_bsf_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("BSF")
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
        str(tmp_path / "Fe*.agr"),
    )


def test_jxc_sp_srel_declares_scalar_and_plot_outputs(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.MODE.MODE = "SP-SREL"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe*Jij*"),
        str(tmp_path / "Fe*.agr"),
    )
