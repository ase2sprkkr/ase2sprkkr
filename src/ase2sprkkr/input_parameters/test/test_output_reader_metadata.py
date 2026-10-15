import pytest

from ..input_parameters import InputParameters
from ..task_definition import TaskDefinition
from ...common.module_registry import RegistryEntryNotFound
from ...outputs.reader_registry import reader_registry
from ...outputs.task_result import KkrOutputReader


def _reader_for_task(task):
    if task.output_reader is not None:
        return reader_registry[task.output_reader]
    return reader_registry.class_for_task(task.name)


def test_all_registered_tasks_have_resolvable_output_reader_metadata():
    for task in InputParameters.task_registry.values():
        assert issubclass(_reader_for_task(task), KkrOutputReader)


def test_explicit_output_reader_name_is_strict():
    task = TaskDefinition(
        "SCF",
        name="BROKEN-READER-MAPPING",
        output_reader="NO-SUCH-READER",
    )

    with pytest.raises(RegistryEntryNotFound):
        reader_registry[task.output_reader]
