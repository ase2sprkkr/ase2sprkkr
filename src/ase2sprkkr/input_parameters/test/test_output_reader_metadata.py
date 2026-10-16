import pytest

from ..input_parameters import InputParameters
from ..task_definition import TaskDefinition
from ...common.module_registry import RegistryEntryNotFound
from ...outputs.reader_registry import reader_registry
from ...outputs.task_result import KkrOutputReader


def test_all_registered_tasks_have_resolvable_output_reader_metadata():
    for task in InputParameters.task_registry.values():
        assert issubclass(reader_registry.class_for_task(task), KkrOutputReader)


def test_explicit_output_reader_name_is_strict():
    task = TaskDefinition(
        "SCF",
        name="BROKEN-READER-MAPPING",
        output_reader="NO-SUCH-READER",
    )

    with pytest.raises(RegistryEntryNotFound):
        reader_registry.class_for_task(task)


def test_bsf_task_names_resolve_explicit_shared_reader():
    assert reader_registry.class_for_task_name("BSFEK").__name__ == "BsfOutputReader"
    assert reader_registry.class_for_task_name("BSFKK").__name__ == "BsfOutputReader"
