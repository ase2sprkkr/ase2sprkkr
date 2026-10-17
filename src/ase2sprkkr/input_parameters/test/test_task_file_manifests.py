import pytest

from ase2sprkkr.common.module_registry import RegistryEntryNotFound
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.task_definition import TaskDefinition, TaskFile
from ase2sprkkr.result_files.definition_registry import result_file_definition_registry


def test_dos_output_is_namespaced_by_dataset(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.dos"),)


def test_dos_output_reuses_result_file_definition_metadata():
    task = InputParameters.task_registry["DOS"]
    file = task.output_files[0]

    assert file.pathname == "{dataset}.dos"
    assert file.type == "DOS"
    assert file.info is None
    definition = result_file_definition_registry.definition_for_task_file(file)
    assert definition is result_file_definition_registry["DOS"]
    assert result_file_definition_registry.info_for_task_file(file) == definition.info


def test_task_file_info_overrides_result_file_definition_info():
    file = TaskFile("result.dos", type="DOS", info="Task-specific meaning")

    assert result_file_definition_registry.info_for_task_file(file) == "Task-specific meaning"


def test_task_file_type_is_a_strict_result_definition_reference():
    file = TaskFile("result.unknown", type="NO-SUCH-TYPE")

    with pytest.raises(RegistryEntryNotFound):
        result_file_definition_registry.definition_for_task_file(file)


def test_scf_declares_new_potential_as_dotted_template(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = "Fe.pot"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.pot_new"),)


def test_missing_dotted_template_value_skips_optional_file(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = None

    assert parameters.output_files(str(tmp_path)) == ()


def test_empty_dotted_template_value_is_rejected(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.POTFIL = ""

    with pytest.raises(ValueError, match=r"CONTROL\.POTFIL.*must not be empty"):
        parameters.output_files(str(tmp_path))


def test_scf_fullpot_shape_function_is_callback_cache(tmp_path):
    parameters = InputParameters.create("SCF")
    parameters.SCF.FULLPOT = True
    parameters.CONTROL.SFNFIL = "Fe.sfn"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.sfn"),)

    parameters.SCF.FULLPOT = False
    assert parameters.cache_files(str(tmp_path)) == ()


def test_custom_task_metadata_can_reuse_registered_schema(tmp_path):
    parameters = InputParameters.create("SCF")
    task = TaskDefinition(
        "SCF",
        name="SCF",
        executable="fixed-kkrscf",
        output_files=(TaskFile("different.dat"),),
    )
    registry = InputParameters.input_parameters_registry

    assert registry.definition_for_task(task) is registry["SCF"]
    assert registry.files_for_task(task, parameters, str(tmp_path), "output") == (
        str(tmp_path / "different.dat"),
    )


def test_callback_programming_error_is_not_hidden(tmp_path):
    parameters = InputParameters.create("SCF")

    def broken(_parameters):
        raise RuntimeError("broken task file callback")

    with pytest.raises(RuntimeError, match="broken task file callback"):
        InputParameters.input_parameters_registry._resolve_declared_files(
            parameters, str(tmp_path), (TaskFile(broken),)
        )


def test_callback_empty_path_is_rejected(tmp_path):
    parameters = InputParameters.create("SCF")

    with pytest.raises(ValueError, match="pathname must not be empty"):
        InputParameters.input_parameters_registry._resolve_declared_files(
            parameters, str(tmp_path), (TaskFile(lambda _parameters: ""),)
        )


def test_xas_reuses_dataset_tau_directory_as_optional_cache(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.cache_files(str(tmp_path)) == (str(tmp_path / "Fe.tau.dir"),)


def test_tau_callback_respects_taupath_and_taudir(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.TAU.TAUPATH = "tau-cache"
    parameters.TAU.TAUDIR = "custom.tau"

    assert parameters.cache_files(str(tmp_path)) == (
        str(tmp_path / "tau-cache" / "custom.tau"),
    )


def test_tau_callback_preserves_absolute_taudir(tmp_path):
    parameters = InputParameters.create("XAS")
    tau = tmp_path / "shared.tau"
    parameters.TAU.TAUDIR = str(tau)
    parameters.TAU.TAUPATH = "ignored-for-absolute-path"

    assert parameters.cache_files(str(tmp_path / "calculation")) == (str(tau),)


def test_empty_taudir_is_rejected(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.TAU.TAUDIR = ""

    with pytest.raises(ValueError, match="pathname must not be empty"):
        parameters.cache_files(str(tmp_path))


def test_empty_taupath_is_rejected(tmp_path):
    parameters = InputParameters.create("XAS")
    parameters.TAU.TAUDIR = "custom.tau"
    parameters.TAU.TAUPATH = ""

    with pytest.raises(ValueError, match="pathname must not be empty"):
        parameters.cache_files(str(tmp_path))


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

    parameters.TASK.SETFLUCT = "MLIN"
    assert parameters.input_files(str(tmp_path)) == ()


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

    tensor, scalar, plot = parameters.task_definition.output_files
    assert tensor.type == "JXC"
    assert scalar.type == "JXC"
    assert plot.type is None


def test_jxc_sp_srel_declares_scalar_and_plot_outputs(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = "Fe"
    parameters.MODE.MODE = "SP-SREL"

    assert parameters.output_files(str(tmp_path)) == (
        str(tmp_path / "Fe*Jij*"),
        str(tmp_path / "Fe*.agr"),
    )


def test_jxc_empty_dataset_is_rejected(tmp_path):
    parameters = InputParameters.create("JXC")
    parameters.CONTROL.DATASET = ""
    parameters.MODE.MODE = "FREL"

    with pytest.raises(ValueError, match="dataset.*must not be empty"):
        parameters.output_files(str(tmp_path))
