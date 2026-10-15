from types import SimpleNamespace

from .. import calculator as calculator_module
from ..calculator import SPRKKR


def test_create_run_forwards_persistent_policy(monkeypatch, tmp_path):
    calculator = SPRKKR(directory=tmp_path)
    persistent = object()
    prepared = object()
    runner = object()
    captured = {}

    monkeypatch.setattr(
        calculator,
        "_prepare_calculation",
        lambda **kwargs: prepared,
    )

    def create_sprkkr_run(value, **kwargs):
        assert value is prepared
        captured.update(kwargs)
        return runner

    monkeypatch.setattr(calculator_module, "_create_sprkkr_run", create_sprkkr_run)

    assert calculator.create_run(
        directory=tmp_path,
        persistent=persistent,
    ) is runner
    assert captured["persistent"] is persistent


def test_calculate_forwards_persistent_policy(monkeypatch):
    calculator = SPRKKR()
    persistent = object()
    captured = {}
    result = SimpleNamespace(task_name="DOS")

    class Runner:
        def run(self):
            return result

    def create_run(**kwargs):
        captured.update(kwargs)
        return Runner()

    monkeypatch.setattr(calculator, "create_run", create_run)

    assert calculator.calculate(persistent=persistent) is result
    assert captured["persistent"] is persistent
