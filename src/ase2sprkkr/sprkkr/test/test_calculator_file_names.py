from types import SimpleNamespace

from ..calculator import SPRKKR


def test_create_run_forwards_explicit_file_names(monkeypatch):
    calculator = SPRKKR(directory=False, empty_spheres=False)
    captured = {}
    runner = object()

    def create_run(**kwargs):
        captured["run"] = kwargs
        return runner

    prepared = SimpleNamespace(create_run=create_run)

    def prepare(**kwargs):
        captured.update(kwargs)
        return prepared

    monkeypatch.setattr(calculator, "prepare", prepare)

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
