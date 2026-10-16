import pytest

from ..input_parameters import InputParameters
from ..task_definition import TaskDefinition
from ...common.module_registry import RegistryEntryNotFound
from ...outputs.reader_registry import reader_registry
from ...outputs.readers.default import DefaultOutputReader
from ...outputs.task_result import KkrOutputReader


def test_all_registered_tasks_have_resolvable_output_reader_metadata():
    for task in InputParameters.task_registry.values():
        assert issubclass(reader_registry.class_for_task(task), KkrOutputReader)


def test_specialized_builtin_readers_are_declared_at_the_right_level():
    expected = {
        "SCF": "ScfOutputReader",
        "DOS": "DosOutputReader",
        "JXC": "JxcOutputReader",
        "BSFEK": "BsfOutputReader",
        "BSFKK": "BsfOutputReader",
        "ARPES": "ArpesOutputReader",
        "AIPES": "SpecOutputReader",
        "BAND": "SpecOutputReader",
        "SPLEED": "SpecOutputReader",
    }

    for task_name, reader_name in expected.items():
        assert reader_registry.class_for_task_name(task_name).__name__ == reader_name


def test_generic_builtin_tasks_use_executable_default_reader():
    for task_name in ("PHAGEN", "TORQUE", "XAS", "GILBERT"):
        assert reader_registry.class_for_task_name(task_name) is DefaultOutputReader


def test_explicit_output_reader_name_is_strict_even_for_custom_executable():
    task = TaskDefinition(
        "SCF",
        name="BROKEN-READER-MAPPING",
        executable="my-custom-sprkkr",
        output_reader="NO-SUCH-READER",
    )

    with pytest.raises(RegistryEntryNotFound):
        reader_registry.class_for_task(task)


def test_explicit_output_reader_does_not_require_registered_executable():
    task = TaskDefinition(
        "SCF",
        name="CUSTOM-SCF",
        executable="my-custom-sprkkr",
        output_reader="SCF",
    )

    assert reader_registry.class_for_task(task).__name__ == "ScfOutputReader"


def test_unknown_executable_without_explicit_reader_is_not_guessed():
    task = TaskDefinition(
        "SCF",
        name="CUSTOM",
        executable="my-custom-sprkkr",
    )

    with pytest.raises(RegistryEntryNotFound, match="my-custom-sprkkr"):
        reader_registry.class_for_task(task)


def test_bsf_task_names_resolve_explicit_shared_reader():
    assert reader_registry.class_for_task_name("BSFEK").__name__ == "BsfOutputReader"
    assert reader_registry.class_for_task_name("BSFKK").__name__ == "BsfOutputReader"
