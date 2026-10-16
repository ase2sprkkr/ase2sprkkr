from types import SimpleNamespace

from ...input_parameters.input_parameters import InputParameters
from ...input_parameters.task_definition import TaskDefinition
from ...sprkkr.executables import executables
from ..reader_registry import OutputReaderRegistry, reader_registry
from ..readers.default import DefaultOutputReader


def test_default_reader_registry_dependencies_are_explicit():
    assert reader_registry._tasks is InputParameters.task_registry
    assert reader_registry._executables is executables


def test_custom_reader_registry_uses_injected_executables():
    registry = OutputReaderRegistry(
        executables={
            "custom-kkr": SimpleNamespace(default_output_reader="DEFAULT")
        }
    )
    task = TaskDefinition(
        "SCF",
        name="CUSTOM",
        executable="custom-kkr",
    )

    assert registry.class_for_task(task) is DefaultOutputReader


def test_custom_reader_registry_uses_injected_task_namespace():
    task = TaskDefinition(
        "SCF",
        name="CUSTOM",
        output_reader_name="SCF",
    )
    registry = OutputReaderRegistry(tasks={"CUSTOM": task})

    assert registry.class_for_task_name("CUSTOM") is registry["SCF"]