from .. import calculator as calculator_module
from ..calculator import SPRKKR


def test_create_run_forwards_explicit_file_names(monkeypatch):
    calculator = SPRKKR(directory=False, empty_spheres=False)
    captured = {}
    prepared = object()
    runner = object()

    def prepare(**kwargs):
        captured.update(kwargs)
        return prepared

    def create_sprkkr_run(value, **kwargs):
        assert value is prepared
        return runner

    monkeypatch.setattr(calculator, "_prepare_calculation", prepare)
    monkeypatch.setattr(calculator_module, "_create_sprkkr_run", create_sprkkr_run)

    assert (
        calculator.create_run(
            input_file="custom.inp",
            potential_file="custom.pot",
            output_file="custom.out",
            directory=False,
        )
        is runner
    )

    assert captured["input_file"] == "custom.inp"
    assert captured["potential_file"] == "custom.pot"
    assert captured["output_file"] == "custom.out"
    assert captured["directory"] is False
