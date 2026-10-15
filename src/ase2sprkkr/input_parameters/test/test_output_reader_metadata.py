import pytest

from ..input_parameters import InputParameters
from ..task_definition import TaskDefinition
from ...common.module_registry import RegistryEntryNotFound
from ...outputs.task_result import KkrOutputReader


def test_all_registered_tasks_resolve_output_reader():
    for task in InputParameters.task_registry.values():
        assert issubclass(task.output_reader_class(), KkrOutputReader)


def test_explicit_output_reader_name_is_strict():
    task = TaskDefinition(
        "SCF",
        name="BROKEN-READER-MAPPING",
        output_reader="NO-SUCH-READER",
    )

    with pytest.raises(RegistryEntryNotFound):
        task.output_reader_class()
