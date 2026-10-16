import pytest

from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..executables import ExecutableDefinition, ExecutableRegistry, executable_registry


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


def test_executable_registry_can_be_extended_independently_of_task_registry():
    registry = ExecutableRegistry()
    definition = ExecutableDefinition("custom", default_output_reader="SPEC")

    assert registry.register(definition) is definition
    assert registry["custom"] is definition


def test_executable_registry_rejects_duplicate_definitions():
    registry = ExecutableRegistry((ExecutableDefinition("custom"),))

    with pytest.raises(RegistryDefinitionError, match="Duplicate executable"):
        registry.register(ExecutableDefinition("custom"))


def test_executable_registry_reports_unknown_executable():
    with pytest.raises(RegistryEntryNotFound, match="unknown"):
        executable_registry["unknown"]


def test_builtin_executable_reader_defaults():
    assert executable_registry["kkrscf"].default_output_reader == "DEFAULT"
    assert executable_registry["kkrgen"].default_output_reader == "DEFAULT"
    assert executable_registry["kkrspec"].default_output_reader == "SPEC"
    assert executable_registry["kkrchi"].default_output_reader == "DEFAULT"
