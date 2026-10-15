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
    assert parameters.requested_task_definition is parameters.task_definition
    assert parameters._definition is parameters.task_definition.input_parameters_definition
    assert parameters.task_definition.executable == "kkrscf"
    assert parameters.task_definition.mpi is True


@pytest.mark.parametrize("copy_values", [False, True])
def test_copy_preserves_requested_task_binding(copy_values):
    parameters = InputParameters.create("BSFEK")

    copied = parameters.copy(copy_values=copy_values)

    assert copied.requested_task_definition is parameters.requested_task_definition
    assert copied.task_definition is parameters.task_definition
    assert copied.task_name == "BSFEK"
    assert copied._definition is copied.task_definition.input_parameters_definition
    np.testing.assert_array_equal(copied.ENERGY.NE(), [200])


def test_shared_schema_distinguishes_requested_and_effective_task():
    parameters = InputParameters.create("BSFEK")
    requested = parameters.requested_task_definition

    parameters.ENERGY.NE = [1]

    assert requested is InputParameters.task_registry["BSFEK"]
    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is InputParameters.task_registry["BSFKK"]
    assert parameters.requested_task_name == "BSFEK"
    assert parameters.task_name == "BSFKK"


def test_change_task_changes_requested_task_and_schema_atomically():
    parameters = InputParameters.create("DOS")

    parameters.change_task("SCF")

    assert parameters.requested_task_definition is InputParameters.task_registry["SCF"]
    assert parameters.task_name == "SCF"
    assert parameters.task_definition.name == "SCF"
    assert parameters._definition is parameters.task_definition.input_parameters_definition


def test_shared_spec_task_initializes_discriminator():
    parameters = InputParameters.create("SPLEED")
    requested = parameters.requested_task_definition

    assert parameters.TASK.TASK() == "SPLEED"
    assert parameters.task_definition is requested

    parameters.TASK.TASK = "AIPES"
    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is InputParameters.task_registry["AIPES"]


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


@pytest.mark.parametrize(
    ("task_name", "reader_name"),
    [
        ("ARPES", "ArpesOutputReader"),
        ("AIPES", "SpecOutputReader"),
        ("SPLEED", "SpecOutputReader"),
        ("BAND", "SpecOutputReader"),
        ("BSFEK", "BsfOutputReader"),
        ("BSFKK", "BsfOutputReader"),
    ],
)
def test_effective_task_owns_reader_contract(task_name, reader_name):
    task = InputParameters.task_registry[task_name]

    assert task.output_reader_class().__name__ == reader_name
