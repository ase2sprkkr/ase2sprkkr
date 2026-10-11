from types import SimpleNamespace

from ..calculator import SPRKKR


def test_run_forwards_explicit_file_names(monkeypatch):
    calculator = SPRKKR(directory=False, empty_spheres=False)
    captured = {}

    class InputParameters:
        TASK = SimpleNamespace(TASK=lambda: "DOS")

        @staticmethod
        def run_process(*args, **kwargs):
            return SimpleNamespace()

    def save_input(**kwargs):
        captured.update(kwargs)
        return InputParameters(), "prepared.inp", None, "prepared.out"

    monkeypatch.setattr(calculator, "save_input", save_input)

    calculator.run(
        input_file="custom.inp",
        potential_file="custom.pot",
        output_file="custom.out",
        directory=False,
    )

    assert captured["input_file"] == "custom.inp"
    assert captured["potential_file"] == "custom.pot"
    assert captured["output_file"] == "custom.out"
