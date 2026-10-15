from types import SimpleNamespace

from ase.build import bulk

from ...sprkkr import calculator as calculator_module
from ...sprkkr.calculator import SPRKKR


class _Runner:
    def __init__(self, result=None):
        self.result = result
        self.ran = False

    def run(self):
        self.ran = True
        return self.result


def test_create_run_only_prepares_runner(tmp_path, monkeypatch):
    calculator = SPRKKR(
        atoms=bulk("Al"),
        directory=tmp_path,
        mpi=False,
        empty_spheres=False,
    )
    runner = _Runner()
    prepared = object()
    seen = {}

    def prepare(**kwargs):
        seen["prepare"] = kwargs
        return prepared

    def create_run(prepared_calculation, **kwargs):
        seen["prepared"] = prepared_calculation
        seen["run"] = kwargs
        return runner

    monkeypatch.setattr(calculator, "_prepare_calculation", prepare)
    monkeypatch.setattr(calculator_module, "_create_sprkkr_run", create_run)

    created = calculator.create_run(
        input_file="custom.inp",
        potential_file="custom.pot",
        create_subdirs=True,
    )

    assert created is runner
    assert not runner.ran
    assert seen["prepared"] is prepared
    assert seen["prepare"]["input_file"] == "custom.inp"
    assert seen["prepare"]["potential_file"] == "custom.pot"
    assert seen["prepare"]["create_subdirs"] is True
    assert seen["run"]["result_context"] is calculator


def test_calculate_runs_and_updates_state(monkeypatch):
    calculator = SPRKKR(atoms=bulk("Al"), mpi=False, empty_spheres=False)
    new_potential = SimpleNamespace(atoms=None)
    result = SimpleNamespace(task_name="SCF", energy=1.25, potential=new_potential)
    runner = _Runner(result)
    callback_state = []

    monkeypatch.setattr(calculator, "create_run", lambda **kwargs: runner)

    returned = calculator.calculate(
        callback=lambda out: callback_state.append(
            (out, calculator.results["energy"], calculator._potential)
        )
    )

    assert returned is result
    assert runner.ran
    assert calculator.results["energy"] == 1.25
    assert calculator._potential is new_potential
    assert callback_state == [(result, 1.25, new_potential)]


def test_calculate_does_not_adopt_potential_when_explicit_atoms_are_used(monkeypatch):
    calculator = SPRKKR(atoms=bulk("Al"), mpi=False, empty_spheres=False)
    old_potential = calculator.potential
    new_potential = SimpleNamespace(atoms=None)
    result = SimpleNamespace(task_name="SCF", energy=1.25, potential=new_potential)
    runner = _Runner(result)

    monkeypatch.setattr(calculator, "create_run", lambda **kwargs: runner)

    calculator.calculate(atoms=bulk("Al"))

    assert calculator._potential is old_potential
