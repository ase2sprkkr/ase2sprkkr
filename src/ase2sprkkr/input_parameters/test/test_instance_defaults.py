from types import SimpleNamespace

import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import (
    InputParametersDefinition,
    InputSectionDefinition,
    InputValueDefinition,
)


def _repeated_schema():
    return InputParametersDefinition(
        "TEST",
        [
            InputSectionDefinition(
                "ITEM",
                [InputValueDefinition("VALUE", 1)],
                is_repeated=True,
                is_optional=True,
            )
        ],
    )


def test_copy_preserves_empty_repeated_items():
    parameters = _repeated_schema().create_object()
    parameters.ITEM.add()

    copied = parameters.copy()

    assert len(parameters.ITEM) == 1
    assert len(copied.ITEM) == 1


def test_copy_preserves_repeated_item_containing_instance_default():
    parameters = _repeated_schema().create_object()
    item = parameters.ITEM.add()
    item.VALUE.set_default(2)

    copied = parameters.copy()

    assert len(copied.ITEM) == 1
    assert copied.ITEM[0].VALUE() == 2
    assert copied.ITEM[0].VALUE.default_value == 2
    assert not copied.ITEM[0].VALUE.is_set()


def test_empty_reader_directory_is_explicit(tmp_path):
    parameters = InputParameters.create("SCF")
    task = parameters.task_definition_required()

    reader = task.output_reader(
        parameters,
        calculator=SimpleNamespace(directory=tmp_path),
        directory="",
    )

    assert reader.directory == ""


def test_false_reader_directory_is_rejected():
    parameters = InputParameters.create("SCF")
    task = parameters.task_definition_required()

    with pytest.raises(TypeError, match="directory must be a path-like object or None"):
        task.output_reader(parameters, directory=False)
