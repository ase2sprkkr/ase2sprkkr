from types import SimpleNamespace

from ..input_parameters import InputParameters
from ..task_definition import TASKS, TaskDefinition


def test_runtime_input_parameters_are_bound_to_task_definition():
    parameters = InputParameters.create("SCF")

    assert isinstance(parameters.task_definition, TaskDefinition)
    assert parameters.task_definition is TASKS.definition("SCF")
    assert parameters._definition is parameters.task_definition.input_parameters_definition
    assert parameters.task_name == "SCF"


def test_copy_preserves_task_binding():
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "copy-test"

    copied = parameters.copy(copy_values=True)

    assert copied.task_definition is parameters.task_definition
    assert copied._definition is parameters._definition
    assert copied.CONTROL.DATASET() == "copy-test"


def test_change_task_switches_task_and_schema_together():
    parameters = InputParameters.create("SCF")

    parameters.change_task("DOS")

    assert parameters.task_definition is TASKS.definition("DOS")
    assert parameters._definition is parameters.task_definition.input_parameters_definition
    assert parameters.task_name == "DOS"


def test_bsf_aliases_share_schema_but_keep_task_identity():
    bsf = TASKS.definition("BSF")
    bsfek = InputParameters.create("BSFEK")
    bsfkk = InputParameters.create("BSFKK")

    assert bsfek.task_definition is TASKS.definition("BSFEK")
    assert bsfkk.task_definition is TASKS.definition("BSFKK")
    assert bsfek.task_definition is not bsf
    assert bsfkk.task_definition is not bsf
    assert bsfek._definition is bsf.input_parameters_definition
    assert bsfkk._definition is bsf.input_parameters_definition
    assert bsfek.task_name == "BSFEK"
    assert bsfkk.task_name == "BSFKK"


def test_task_policy_is_not_on_input_schema():
    task = TASKS.definition("PHAGEN")
    schema = task.input_parameters_definition

    assert task.executable == "kkrscf"
    assert task.mpi is False
    for attribute in (
        "executable",
        "mpi",
        "result_reader",
        "input_files",
        "output_files",
        "cache_files",
        "input_file_templates",
        "output_file_templates",
        "cache_file_templates",
    ):
        assert not hasattr(schema, attribute)


def test_dos_output_manifest_belongs_to_task(tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "Fe"

    assert parameters.output_files(str(tmp_path)) == (str(tmp_path / "Fe.dos"),)


def test_jxc_continuation_policy_belongs_to_target_task():
    parameters = InputParameters.create("JXC")
    result = SimpleNamespace(
        last_iteration=SimpleNamespace(
            energy=SimpleNamespace(EMIN=lambda: -0.42),
        )
    )

    parameters.task_definition.prepare_parameters_from_result(parameters, result)

    assert parameters.ENERGY.EMIN() == -0.42
