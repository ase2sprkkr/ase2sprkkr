from pathlib import Path

from ase.build import bulk

from ..calculator import SPRKKR


def _calculator(tmp_path):
    return SPRKKR(
        atoms=bulk("Li"),
        directory=tmp_path,
        output_file="default.out",
        empty_spheres=False,
        mpi=False,
    )


def test_save_input_without_output_has_stable_four_item_result(tmp_path):
    calculator = _calculator(tmp_path)

    input_parameters, input_file, potential_file, output_file = calculator.save_input(
        output_file=False
    )

    assert input_parameters is calculator.input_parameters
    assert Path(input_file).is_file()
    assert Path(potential_file).is_file()
    assert output_file is None
    assert not (tmp_path / "default.out").exists()


def test_create_run_preserves_explicitly_disabled_output(tmp_path, monkeypatch):
    calculator = _calculator(tmp_path)
    captured = {}

    def fake_run_process(self, calculator, input_file, output_file, **kwargs):
        captured["output_file"] = output_file
        input_file.close()
        return object()

    monkeypatch.setattr(type(calculator.input_parameters), "run_process", fake_run_process)

    calculator.create_run(output_file=False)

    assert captured["output_file"] is None
    assert not (tmp_path / "default.out").exists()
