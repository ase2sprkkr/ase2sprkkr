from dataclasses import FrozenInstanceError

import pytest

from ase2sprkkr import executables as public_executables
from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..executables import ExecutableDefinition, _ExecutableRegistry, executables


def test_executable_definition_normalizes_reader_name_only():
    definition = ExecutableDefinition(
        "CustomExecutable",
        default_output_reader=" spec ",
        description="Custom executable",
    )

    assert definition.name == "CustomExecutable"
    assert definition.default_output_reader == "SPEC"
    assert definition.description == "Custom executable"


def test_executable_definition_keeps_executable_name_exact():
    definition = ExecutableDefinition(" custom ")

    assert definition.name == " custom "


@pytest.mark.parametrize(
    ("kwargs", "error"),
    (
        ({"name": None}, TypeError),
        ({"name": "   "}, ValueError),
        ({"name": "custom", "default_output_reader": None}, TypeError),
        ({"name": "custom", "default_output_reader": "   "}, ValueError),
        ({"name": "custom", "description": None}, TypeError),
    ),
)
def test_executable_definition_rejects_invalid_metadata(kwargs, error):
    with pytest.raises(error):
        ExecutableDefinition(**kwargs)


def test_executable_registry_accepts_custom_definitions_at_construction():
    definition = ExecutableDefinition("custom", default_output_reader="SPEC")
    registry = _ExecutableRegistry((definition,))

    assert registry["custom"] is definition
    assert not hasattr(registry, "register")


def test_executable_registry_rejects_non_definition_items():
    with pytest.raises(TypeError, match="ExecutableDefinition"):
        _ExecutableRegistry(("custom",))


def test_executable_registry_rejects_duplicate_definitions():
    with pytest.raises(RegistryDefinitionError, match="Duplicate executable"):
        _ExecutableRegistry(
            (ExecutableDefinition("custom"), ExecutableDefinition("custom"))
        )


def test_executable_registry_reports_unknown_executable():
    with pytest.raises(RegistryEntryNotFound, match="unknown"):
        executables["unknown"]


def test_builtin_executable_reader_defaults():
    assert executables["kkrscf"].default_output_reader == "DEFAULT"
    assert executables["kkrgen"].default_output_reader == "DEFAULT"
    assert executables["kkrspec"].default_output_reader == "SPEC"
    assert executables["kkrchi"].default_output_reader == "DEFAULT"


def test_executables_is_public_read_only_metadata():
    assert public_executables is executables
    assert not hasattr(public_executables, "register")

    with pytest.raises(TypeError):
        public_executables["custom"] = ExecutableDefinition("custom")

    with pytest.raises(FrozenInstanceError):
        public_executables["kkrscf"].description = "Changed"
