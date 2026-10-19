"""The removed persistent flag cannot silently reach input parameters."""

from ase2sprkkr.sprkkr.calculator import SPRKKR


def test_project_is_the_supported_execution_scope():
    calculator = SPRKKR()
    assert calculator.project is None
    assert not hasattr(calculator, "persistent")


def test_project_argument_is_not_an_input_parameter():
    project = object()
    assert SPRKKR(project=project).project is project
