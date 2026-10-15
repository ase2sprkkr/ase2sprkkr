import numpy as np

from ..input_parameters import InputParameters
from ...outputs.readers.bsf import BsfOutputReader


def test_requested_task_is_task_definition_object():
    parameters = InputParameters.create("SCF")
    task = InputParameters.task_registry["SCF"]

    assert parameters.requested_task_definition is task
    assert parameters.task_definition is task
    assert parameters.requested_task_name == "SCF"
    assert parameters.task_name == "SCF"


def test_spec_effective_task_changes_but_requested_task_does_not():
    parameters = InputParameters.create("SPLEED")
    requested = InputParameters.task_registry["SPLEED"]

    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is requested
    assert parameters.TASK.TASK() == "SPLEED"

    parameters.TASK.TASK = "ARPES"

    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is InputParameters.task_registry["ARPES"]
    assert parameters.requested_task_name == "SPLEED"
    assert parameters.task_name == "ARPES"


def test_bsf_effective_task_changes_but_requested_task_does_not():
    parameters = InputParameters.create("BSFEK")
    requested = InputParameters.task_registry["BSFEK"]

    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is requested
    np.testing.assert_array_equal(parameters.ENERGY.NE(), [200])

    parameters.ENERGY.NE = [1]

    assert parameters.requested_task_definition is requested
    assert parameters.task_definition is InputParameters.task_registry["BSFKK"]
    assert parameters.requested_task_name == "BSFEK"
    assert parameters.task_name == "BSFKK"


def test_copy_preserves_requested_provenance_and_resolves_current_values():
    parameters = InputParameters.create("BSFEK")
    parameters.ENERGY.NE = [1]

    copied = parameters.copy(copy_values=True)

    assert copied.requested_task_definition is parameters.requested_task_definition
    assert copied.task_definition is InputParameters.task_registry["BSFKK"]
    np.testing.assert_array_equal(copied.ENERGY.NE(), [1])


def test_parsed_shared_schema_has_no_requested_task():
    parameters = InputParameters.create("SPLEED")
    parameters.CONTROL.POTFIL = "test.pot"
    text = parameters.to_string(validate=False)

    parsed = InputParameters.definition("SPLEED").read_from_string(
        text,
        allow_dangerous=True,
    )

    assert parsed.requested_task_definition is None
    assert parsed.task_definition is InputParameters.task_registry["SPLEED"]


def test_change_task_sets_requested_task_and_initializes_discriminator():
    parameters = InputParameters.create("DOS")

    parameters.change_task("SPLEED")

    assert parameters.requested_task_definition is InputParameters.task_registry["SPLEED"]
    assert parameters.task_definition is parameters.requested_task_definition
    assert parameters.TASK.TASK() == "SPLEED"


def test_change_task_initializes_bsf_mode_from_requested_task():
    parameters = InputParameters.create("BSFKK")

    parameters.change_task("BSFEK")

    assert parameters.requested_task_definition is InputParameters.task_registry["BSFEK"]
    assert parameters.task_definition is parameters.requested_task_definition
    np.testing.assert_array_equal(parameters.ENERGY.NE(), [200])


def test_bsf_tasks_share_schema_but_keep_task_specific_reader_contract():
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek.input_parameters_definition is bsfkk.input_parameters_definition
    assert bsfek.output_reader_class() is BsfOutputReader
    assert bsfkk.output_reader_class() is BsfOutputReader
