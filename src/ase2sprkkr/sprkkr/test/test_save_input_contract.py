from pathlib import Path

from ase.build import bulk
import pytest

from ...input_parameters.input_parameters import InputParameters
from ...potentials.potentials import Potential
from ..calculation import PreparedCalculation
from ..calculator import SPRKKR


@pytest.fixture(scope="module")
def li_bulk():
    # These tests need independent calculator state, not six identical
    # constructions of the same pristine Li lattice.
    return bulk("Li")


def _calculator(tmp_path, atoms, *, output_file="default.out"):
    return SPRKKR(
        atoms=atoms.copy(),
        directory=tmp_path,
        output_file=output_file,
        empty_spheres=False,
        mpi=False,
    )


def test_save_input_without_output_has_stable_four_item_result(tmp_path, li_bulk):
    calculator = _calculator(tmp_path, li_bulk)
    defaults = calculator.input_parameters
    original_potfil = defaults.CONTROL.POTFIL.result

    input_parameters, input_file, potential_file, output_file = calculator.save_input(
        output_file=False
    )

    assert input_parameters is not defaults
    assert defaults.CONTROL.POTFIL.result == original_potfil
    assert Path(input_file).is_file()
    assert Path(potential_file).is_file()
    assert output_file is None
    assert not (tmp_path / "default.out").exists()


@pytest.fixture
def lightweight_serialization(monkeypatch):
    """Preserve real file setup, avoiding unrelated POT/IPD serialization."""
    calls = []

    def save_potential(self, file, *args, **kwargs):
        calls.append("potential")
        file.write("# potential routing test\n")

    def save_parameters(self, file, *args, **kwargs):
        calls.append("input")
        file.write("# input routing test\n")

    monkeypatch.setattr(Potential, "save_to_file", save_potential)
    monkeypatch.setattr(InputParameters, "save_to_file", save_parameters)
    yield
    assert sorted(calls) == ["input", "potential"]


def test_create_run_preserves_explicitly_disabled_output(
    tmp_path, monkeypatch, li_bulk, lightweight_serialization
):
    calculator = _calculator(tmp_path, li_bulk)
    captured = {}
    runner = object()

    def create_run(prepared, **kwargs):
        captured["output_file"] = prepared.output_file
        prepared.close()
        return runner

    monkeypatch.setattr(PreparedCalculation, "create_run", create_run)

    assert calculator.create_run(output_file=False) is runner

    assert captured["output_file"] is None
    assert not (tmp_path / "default.out").exists()


def test_create_run_respects_calculator_level_disabled_output(
    tmp_path, monkeypatch, li_bulk, lightweight_serialization
):
    calculator = _calculator(tmp_path, li_bulk, output_file=False)
    captured = {}
    runner = object()

    def create_run(prepared, **kwargs):
        captured["output_file"] = prepared.output_file
        prepared.close()
        return runner

    monkeypatch.setattr(PreparedCalculation, "create_run", create_run)

    assert calculator.create_run() is runner
    assert captured["output_file"] is None


def test_relative_input_parameter_file_is_resolved_against_work_directory(tmp_path, li_bulk):
    source = _calculator(tmp_path, li_bulk)
    _, input_name, potential_name, _ = source.save_input(
        input_file="source.inp",
        potential_file="source.pot",
        output_file=False,
    )

    assert Path(input_name) == tmp_path / "source.inp"
    assert Path(potential_name) == tmp_path / "source.pot"

    calculator = SPRKKR(
        directory=tmp_path,
        potential=False,
        output_file=False,
        empty_spheres=False,
        mpi=False,
    )
    parameters, reused_input, reused_potential, output = calculator.save_input(
        input_parameters="source.inp",
        potential=False,
        output_file=False,
    )

    assert Path(reused_input) == tmp_path / "source.inp"
    assert Path(reused_potential) == tmp_path / parameters.CONTROL.POTFIL.result
    assert Path(reused_potential).is_file()
    assert output is None


def test_save_input_rejects_temporary_work_directory(tmp_path, li_bulk):
    calculator = _calculator(tmp_path, li_bulk)

    with pytest.raises(ValueError, match="persistent working directory"):
        calculator.save_input(directory=False)


def test_return_files_requires_persistent_potential_path(tmp_path, li_bulk):
    calculator = _calculator(tmp_path, li_bulk)

    with pytest.raises(ValueError, match="temporary POTFIL"):
        calculator.save_input(
            input_file=False,
            potential_file=False,
            output_file=False,
            return_files=True,
        )
