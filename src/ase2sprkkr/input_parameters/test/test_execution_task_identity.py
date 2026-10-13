from types import SimpleNamespace

from ..input_parameters import InputParameters


def test_execution_task_name_uses_task_from_concrete_input():
    parameters = SimpleNamespace(
        task_name="DEFINITION",
        TASK=SimpleNamespace(TASK=lambda: "ACTUAL"),
    )

    assert InputParameters._execution_task_name.fget(parameters) == "ACTUAL"


def test_execution_task_name_falls_back_when_task_is_not_encoded():
    parameters = SimpleNamespace(task_name="DEFINITION")

    assert InputParameters._execution_task_name.fget(parameters) == "DEFINITION"
