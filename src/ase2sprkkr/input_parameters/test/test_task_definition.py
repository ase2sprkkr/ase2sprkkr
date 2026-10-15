import numpy as np
import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


def _schema_for(task):
    return InputParameters.task_registry.input_parameters_definition(task)


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
    assert parameters._definition is _schema_for(parameters.task_definition)
    assert parameters.task_definition.executable == "kkrscf"
    assert parameters.task_definition.mpi is True


@pytest.mark.parametrize("copy_values", [False, True])
def test_copy_preserves_task_binding(copy_values):
    parameters = InputParameters.create("BSFEK")

    copied = parameters.copy(copy_values=copy_values)

    assert copied.task_definition is parameters.task_definition
    assert copied.task_name == "BSFEK"
    assert copied._definition is _schema_for(copied.task_definition)
    np.testing.assert_array_equal(copied.ENERGY.NE(), [200])


def test_change_task_changes_task_and_schema_atomically():
    parameters = InputParameters.create("DOS")

    parameters.change_task("SCF")

    assert parameters.task_name == "SCF"
    assert parameters.task_definition.name == "SCF"
    assert parameters._definition is _schema_for(parameters.task_definition)


def test_direct_schema_objects_are_detached():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})
    parameters = schema.create_object()

    assert parameters.task_definition is None
    assert parameters.task_name == schema.name
    with pytest.raises(RuntimeError, match="not associated"):
        parameters.is_mpi()


def test_task_definition_requires_schema_name():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})

    with pytest.raises(TypeError, match="must be a string"):
        TaskDefinition(schema)


def test_effective_task_follows_shared_schema_values():
    parameters = InputParameters.create("ARPES")

    assert parameters.task_name == "ARPES"
    assert parameters.task_definition.name == "ARPES"

    parameters.TASK.TASK = "SPLEED"

    assert parameters.task_name == "ARPES"
    assert parameters.task_definition.name == "SPLEED"


def test_bsf_tasks_declare_shared_output_reader_identifier():
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek.output_reader == "BSF"
    assert bsfkk.output_reader == "BSF"
