from types import SimpleNamespace

import numpy as np
import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ...outputs.reader_registry import reader_registry
from ...outputs.readers.default import DefaultOutputReader
from ...outputs.task_result import KkrOutputReader


def test_task_registry_aliases_share_schema_not_identity():
    bsf = InputParameters.task_registry["BSF"]
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek is not bsf
    assert bsfkk is not bsf
    assert bsfek.input_parameters_definition is bsf.input_parameters_definition
    assert bsfkk.input_parameters_definition is bsf.input_parameters_definition
    assert bsfek.name == "BSFEK"
    assert bsfkk.name == "BSFKK"


def test_task_binding_is_authoritative():
    parameters = InputParameters.create("SCF")

    assert parameters.task_name == "SCF"
    assert parameters._definition is parameters.task_definition.input_parameters_definition
    assert parameters.task_definition.executable == "kkrscf"
    assert parameters.task_definition.mpi is True


@pytest.mark.parametrize("copy_values", [False, True])
def test_copy_preserves_task_binding(copy_values):
    parameters = InputParameters.create("BSFEK")

    copied = parameters.copy(copy_values=copy_values)

    assert copied.task_definition is parameters.task_definition
    assert copied.task_name == "BSFEK"
    assert copied._definition is copied.task_definition.input_parameters_definition
    np.testing.assert_array_equal(copied.ENERGY.NE(), [200])


def test_change_task_changes_task_and_schema_atomically():
    parameters = InputParameters.create("DOS")

    parameters.change_task("SCF")

    assert parameters.task_name == "SCF"
    assert parameters.task_definition.name == "SCF"
    assert parameters._definition is parameters.task_definition.input_parameters_definition


def test_direct_schema_objects_are_detached():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})
    parameters = schema.create_object()

    assert parameters.task_definition is None
    assert parameters.task_name == schema.name
    with pytest.raises(RuntimeError, match="not bound"):
        parameters.is_mpi()


def test_jxc_continuation_is_task_behavior():
    parameters = InputParameters.create("JXC")
    result = SimpleNamespace(
        last_iteration=SimpleNamespace(
            energy=SimpleNamespace(EMIN=lambda: -0.35)
        )
    )

    parameters.task_definition.prepare_from_result(parameters, result)

    assert parameters.ENERGY.EMIN() == -0.35


def test_reader_discovery_uses_registry_and_default_fallback():
    assert reader_registry["SCF"].__name__ == "ScfOutputReader"
    assert KkrOutputReader.class_for_task("unknown") is DefaultOutputReader


def test_task_reader_respects_concrete_execution_task():
    parameters = InputParameters.create("ARPES")
    task = parameters.task_definition
    assert task._execution_task_name(parameters) == "ARPES"
    assert task.output_reader_class(parameters).__name__ == "ArpesOutputReader"

    parameters.TASK.TASK = "SPLEED"
    assert task._execution_task_name(parameters) == "SPLEED"
    assert task.output_reader_class(parameters) is DefaultOutputReader
    assert task.output_reader_class(execution_task="SPLEED") is DefaultOutputReader
