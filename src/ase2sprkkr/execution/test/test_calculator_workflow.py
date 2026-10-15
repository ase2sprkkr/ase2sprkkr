from types import SimpleNamespace

from ase.build import bulk

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
    seen = {}

    class Prepared:
        def create_run(self, **kwargs):
            seen["run"] = kwargs
            return runner

    prepared = Prepared()

    def prepare(**kwargs):
        seen["prepare"] = kwargs
        return prepared

    monkeypatch.setattr(calculator, "prepare", prepare)

    created = calculator.create_run(
        input_file="custom.inp",
        potential_file="custom.pot",
        create_subdirs=True,
    )

    assert created is runner
    assert not runner.ran
    assert seen["prepare"]["input_file"] == "custom.inp"
    assert seen["prepare"]["potential_file"] == "custom.pot"
    assert seen["prepare"]["create_subdirs"] is True
    assert seen["run"]["result_context"] is calculator
    assert seen["run"]["create_subdirs"] is True


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
