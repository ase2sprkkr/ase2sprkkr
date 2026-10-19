from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, TaskFile
from ...outputs.readers.default import DefaultOutputReader
from ...sprkkr.configuration import ConfigurationFile


def _schema_for(task):
    return InputParameters.input_parameters_registry.definition_for_task(task)


def test_bsf_task_views_share_schema_not_identity():
    schema = InputParameters.input_parameters_registry["BSF"]
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek is not bsfkk
    assert _schema_for(bsfek) is schema
    assert _schema_for(bsfkk) is schema
    assert bsfek.name == "BSFEK"
    assert bsfkk.name == "BSFKK"


def test_task_binding_is_authoritative():
    parameters = InputParameters.create("SCF")

    assert parameters.task_name == "SCF"
    assert parameters.definition_name == "SCF"
    assert parameters._definition is _schema_for(parameters.task_definition)
    assert parameters.task_definition.executable == "kkrscf"
    assert parameters.task_definition.mpi is True
    assert parameters.task_definition_required() is parameters.task_definition


@pytest.mark.parametrize("copy_values", [False, True])
def test_copy_preserves_task_binding(copy_values):
    parameters = InputParameters.create("BSFEK")

    copied = parameters.copy(copy_values=copy_values)

    assert copied.task_definition is parameters.task_definition
    assert copied.task_name == "BSFEK"
    assert copied.requested_task_name == "BSFEK"
    assert copied.definition_name == "BSF"
    assert copied._definition is _schema_for(copied.task_definition)
    np.testing.assert_array_equal(copied.ENERGY.NE(), [200])
    assert not copied.ENERGY.NE.is_set()


def test_change_task_changes_task_and_schema_atomically():
    parameters = InputParameters.create("DOS")

    parameters.change_task("SCF")

    assert parameters.task_name == "SCF"
    assert parameters.requested_task_name == "SCF"
    assert parameters.definition_name == "SCF"
    assert parameters.task_definition.name == "SCF"
    assert parameters._definition is _schema_for(parameters.task_definition)


def test_change_task_retargets_shared_schema_default():
    parameters = InputParameters.create("SPLEED")

    assert parameters.TASK.TASK() == "SPLEED"
    assert not parameters.TASK.TASK.is_set()

    parameters.change_task("AIPES")

    assert parameters.requested_task_name == "AIPES"
    assert parameters.task_name == "AIPES"
    assert parameters.TASK.TASK() == "AIPES"
    assert not parameters.TASK.TASK.is_set()
    assert parameters.task_definition.name == "AIPES"


def test_change_task_retain_values_preserves_explicit_shared_task_value():
    parameters = InputParameters.create("ARPES")
    parameters.TASK.TASK = "SPLEED"

    parameters.change_task("AIPES", retain_values=True)

    assert parameters.requested_task_name == "AIPES"
    assert parameters.task_name == "SPLEED"
    assert parameters.TASK.TASK.is_set()
    assert parameters.TASK.TASK() == "SPLEED"
    assert parameters.task_definition.name == "SPLEED"
    assert parameters._definition is _schema_for(parameters.task_definition)


def test_change_task_retain_values_preserves_explicit_old_default():
    parameters = InputParameters.create("ARPES")
    parameters.TASK.TASK = "ARPES"
    assert parameters.TASK.TASK.is_set()

    parameters.change_task("SPLEED", retain_values=True)

    assert parameters.requested_task_name == "SPLEED"
    assert parameters.task_name == "ARPES"
    assert parameters.TASK.TASK.is_set()
    assert parameters.TASK.TASK() == "ARPES"
    assert parameters._definition is _schema_for(parameters.task_definition)


def test_direct_schema_objects_are_detached():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})
    parameters = schema.create_object()

    assert parameters.task_definition is None
    assert parameters.task_name is None
    assert parameters.requested_task_name is None
    with pytest.raises(RuntimeError, match="not associated with an SPR-KKR task"):
        parameters.task_definition_required()
    assert parameters.copy().task_name is None
    assert parameters.copy().requested_task_name is None


def test_canonical_schema_identity_does_not_imply_requested_task():
    schema = InputParameters.input_parameters_registry["SCF"]
    parameters = schema.create_object()

    assert parameters.task_definition is InputParameters.task_registry["SCF"]
    assert parameters.task_name == "SCF"
    assert parameters.definition_name == "SCF"
    assert parameters.requested_task_name is None
    assert not hasattr(parameters, "execution_task_name")


def test_task_definition_requires_schema_name():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})

    with pytest.raises(TypeError, match="must be a string"):
        TaskDefinition(schema)
    with pytest.raises(ValueError, match="must not be empty"):
        TaskDefinition("   ")


def test_task_definition_output_reader_is_an_identifier():
    with pytest.raises(TypeError, match="output_reader must be a string or None"):
        TaskDefinition("SCF", output_reader=DefaultOutputReader)

    task = TaskDefinition("SCF", output_reader=" bsf ")
    assert task.output_reader_name == "BSF"


def test_task_definition_reader_alias_and_canonical_field_are_unambiguous():
    with pytest.raises(TypeError, match="Specify only one"):
        TaskDefinition(
            "SCF",
            output_reader_name="SCF",
            output_reader="SCF",
        )


def test_task_definition_dataclass_replace_uses_canonical_reader_field():
    task = TaskDefinition("SCF", output_reader="SCF")

    replaced = replace(task, output_reader_name="DOS")

    assert replaced.output_reader_name == "DOS"
    assert replace(task) == task


def test_task_definition_normalizes_identifiers_but_preserves_executable_name():
    task = TaskDefinition(
        " scf ",
        name=" scf ",
        executable=" KkrScf ",
        output_reader=" scf ",
    )

    assert task.input_parameters_definition_name == "SCF"
    assert task.name == "SCF"
    assert task.executable == " KkrScf "
    assert task.output_reader_name == "SCF"


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"name": 1}, TypeError, "name must be a string"),
        ({"name": "   "}, ValueError, "name must not be empty"),
        ({"executable": 1}, TypeError, "executable must be a string"),
        ({"executable": "   "}, ValueError, "executable must not be empty"),
        ({"mpi": "banana"}, TypeError, "mpi must be a bool"),
        ({"mpi": 1}, TypeError, "mpi must be a bool"),
        ({"output_reader": "   "}, ValueError, "output_reader must not be empty"),
        ({"output_reader_name": "   "}, ValueError, "output_reader must not be empty"),
    ],
)
def test_task_definition_rejects_invalid_scalar_metadata(kwargs, error, message):
    with pytest.raises(error, match=message):
        TaskDefinition("SCF", **kwargs)


@pytest.mark.parametrize(
    "field",
    ["input_file_templates", "output_file_templates", "cache_file_templates"],
)
def test_task_definition_materializes_pathlike_templates_as_strings(field):
    task = TaskDefinition("SCF", **{field: (Path("{dataset}.dat"),)})

    assert getattr(task, field) == ("{dataset}.dat",)
    hash(task)


@pytest.mark.parametrize(
    "field",
    ["input_file_templates", "output_file_templates", "cache_file_templates"],
)
def test_task_definition_rejects_non_pathlike_template_elements(field):
    with pytest.raises(TypeError, match=rf"{field}\[0\]"):
        TaskDefinition("SCF", **{field: ([1, 2],)})


@pytest.mark.parametrize(
    "field",
    ["input_file_templates", "output_file_templates", "cache_file_templates"],
)
def test_task_definition_rejects_single_template_instead_of_iterable(field):
    with pytest.raises(TypeError, match="iterable of path-like strings"):
        TaskDefinition("SCF", **{field: "{dataset}.dat"})


def test_task_definition_is_an_immutable_value_object():
    first = TaskDefinition(
        "scf",
        name="scf",
        executable="kkrscf",
        output_reader="scf",
        output_file_templates=("{dataset}.out",),
    )
    second = TaskDefinition(
        "SCF",
        name="SCF",
        executable="kkrscf",
        output_reader_name="SCF",
        output_file_templates=("{dataset}.out",),
    )

    assert first == second
    assert hash(first) == hash(second)
    with pytest.raises(FrozenInstanceError):
        first.output_reader_name = "BSF"


def test_builtin_task_definitions_are_metadata_objects():
    assert all(
        type(task) is TaskDefinition
        for task in InputParameters.task_registry.values()
    )
    for name in (
        "output_reader_class",
        "prepare_from_result",
        "execution_task_name",
    ):
        assert not hasattr(TaskDefinition, name)


def test_input_parameters_use_generic_save_contract():
    assert not hasattr(InputParametersDefinition, "save_hook")
    assert InputParameters.save_to_file is ConfigurationFile.save_to_file


def test_task_name_follows_shared_schema_values():
    parameters = InputParameters.create("ARPES")

    assert parameters.requested_task_name == "ARPES"
    assert parameters.task_name == "ARPES"
    assert parameters.task_definition.name == "ARPES"

    parameters.TASK.TASK = "SPLEED"

    assert parameters.requested_task_name == "ARPES"
    assert parameters.task_name == "SPLEED"
    assert parameters.task_definition.name == "SPLEED"


def test_bsf_tasks_declare_shared_output_reader_identifier():
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek.output_reader_name == "BSF"
    assert bsfkk.output_reader_name == "BSF"


def test_optional_file_template_retains_metadata(tmp_path):
    task = TaskDefinition(
        "SCF",
        input_file_templates=(TaskFile("aux.dat", optional=True),),
        output_file_templates=(TaskFile("{dataset}_optional-*.dat", optional=True),),
    )
    parameters = InputParameters.create("SCF")
    parameters.CONTROL.DATASET = "Fe"
    registry = InputParameters.input_parameters_registry
    for kind in ("input", "output"):
        resolved = registry._resolve_file_templates(
            parameters, str(tmp_path), getattr(task, f"{kind}_file_templates")
        )
        assert len(resolved) == 1
        assert resolved[0].optional is True
    assert task.input_file_templates[0].optional is True


def test_optional_file_metadata_affects_task_identity():
    required = TaskDefinition("SCF", output_file_templates=(TaskFile("x", optional=False),))
    optional = TaskDefinition("SCF", output_file_templates=(TaskFile("x", optional=True),))
    assert required != optional
    assert hash(required) != hash(optional)
    assert optional.output_file_templates[0] != "x"
    assert required.output_file_templates[0].pathname == "x"
    with pytest.raises(AttributeError):
        optional.output_file_templates[0].optional = False


def test_task_execution_rejects_cross_schema_selector_without_retargeting():
    parameters = InputParameters.create("SPLEED")
    parameters.TASK.TASK = "ARPES"
    assert parameters.task_name == "ARPES"
    with pytest.raises(RuntimeError, match="change_task"):
        parameters.task_definition_required()

    parameters.change_task("ARPES")
    assert parameters._definition is _schema_for(parameters.task_definition)
    assert parameters.task_definition_required().name == "ARPES"

def test_optional_file_declaration_pickle_roundtrip():
    import pickle

    task = TaskDefinition("SCF", output_file_templates=(TaskFile("extra.dat", optional=True),))
    restored = pickle.loads(pickle.dumps(task))
    assert restored == task
    assert restored.output_file_templates[0].optional is True
