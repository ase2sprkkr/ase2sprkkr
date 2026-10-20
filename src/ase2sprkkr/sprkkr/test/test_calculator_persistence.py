"""Calculator forwards explicit project scope into prepared execution."""

from types import SimpleNamespace
from ..calculator import SPRKKR


def test_calculate_forwards_project(monkeypatch):
    calculator = SPRKKR()
    project = object()
    captured = {}
    result = SimpleNamespace(task_name="DOS")
    class Runner:
        def run(self):
            return result
    monkeypatch.setattr(calculator, "create_run",
                        lambda **kwargs: captured.update(kwargs) or Runner())
    assert calculator.calculate(project=project) is result
    assert captured["project"] is project
