from types import SimpleNamespace

import numpy as np
import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition
from ...outputs.reader_registry import reader_registry
from ...outputs.readers.default import DefaultOutputReader
from ...outputs.task_result import KkrOutputReader


def test_bsf_task_views_share_schema_not_identity():
    schema = InputParameters.input_parameters_registry["BSF"]
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek is not bsfkk
    assert bsfek.input_parameters_definition is schema
    assert bsfkk.input_parameters_definition is schema
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
    with pytest.raises(RuntimeError, match="not associated"):
        parameters.is_mpi()


def test_task_definition_requires_schema_name():
    schema = InputParametersDefinition.definition_from_dict({"TEST": []})

    with pytest.raises(TypeError, match="must be a string"):
        TaskDefinition(schema)


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


def test_task_reader_is_resolved_from_task_definition_name():
    parameters = InputParameters.create("ARPES")
    requested_task = parameters.task_definition
    assert requested_task.output_reader_class().__name__ == "ArpesOutputReader"

    parameters.TASK.TASK = "SPLEED"
    assert parameters.task_name == "ARPES"
    assert parameters.task_definition.name == "SPLEED"

    assert requested_task.output_reader_class().__name__ == "ArpesOutputReader"
    assert parameters.task_definition.output_reader_class() is DefaultOutputReader


def test_bsf_tasks_explicitly_share_bsf_output_reader():
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek.output_reader == "BSF"
    assert bsfkk.output_reader == "BSF"
    assert bsfek.output_reader_class().__name__ == "BsfOutputReader"
    assert bsfkk.output_reader_class().__name__ == "BsfOutputReader"
