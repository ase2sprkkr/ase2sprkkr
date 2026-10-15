from pathlib import Path
from types import SimpleNamespace

import pytest

from ...potentials.potentials import Potential
from .. import calculator as calculator_module
from ..calculator import SPRKKR


EXAMPLES = Path(__file__).parents[2] / "potentials" / "examples"


def _save_status(tmp_path, *, input_parameters, status, task=None):
    potential = Potential.from_file(EXAMPLES / "FePt.new.pot")
    potential.SCF_INFO.SCFSTATUS = status
    calculator = SPRKKR()

    parameters, _, potential_file, _ = calculator.save_input(
        potential=potential,
        input_parameters=input_parameters,
        task=task,
        directory=tmp_path,
        input_file="task.inp",
        potential_file="task.pot",
        output_file=False,
        empty_spheres=False,
    )

    saved = Potential.from_file(potential_file)
    return (
        parameters.task_name.upper(),
        potential.SCF_INFO.SCFSTATUS(),
        saved.SCF_INFO.SCFSTATUS(),
    )


@pytest.mark.parametrize(
    ("input_parameters", "task", "status", "expected_task", "expected_status"),
    [
        ("SCF", None, "START", "SCF", "START"),
        ("DOS", None, "START", "DOS", "CONVERGED"),
        ("SCF", "DOS", "START", "DOS", "CONVERGED"),
        ("DOS", None, "ITR-BULK", "DOS", "ITR-BULK"),
    ],
)
def test_scf_status_serialization_follows_effective_task(
    tmp_path,
    input_parameters,
    task,
    status,
    expected_task,
    expected_status,
):
    actual_task, in_memory_status, saved_status = _save_status(
        tmp_path,
        input_parameters=input_parameters,
        task=task,
        status=status,
    )

    assert actual_task == expected_task
    assert in_memory_status == expected_status
    assert saved_status == expected_status


def test_calculate_serializes_status_before_mocked_execution(tmp_path, monkeypatch):
    potential = Potential.from_file(EXAMPLES / "FePt.new.pot")
    potential.SCF_INFO.SCFSTATUS = "START"
    result = SimpleNamespace(task_name="DOS")

    class Runner:
        ran = False

        def run(self):
            self.ran = True
            return result

    runner = Runner()

    def fake_create_run(prepared, **kwargs):
        assert Potential.from_file(prepared.potential_file).SCF_INFO.SCFSTATUS() == "CONVERGED"
        prepared.close()
        return runner

    monkeypatch.setattr(calculator_module, "_create_sprkkr_run", fake_create_run)

    calculator = SPRKKR(directory=tmp_path, mpi=False, empty_spheres=False)
    returned = calculator.calculate(
        potential=potential,
        input_parameters="SCF",
        task="DOS",
        input_file="task.inp",
        potential_file="task.pot",
        output_file=False,
    )

    assert returned is result
    assert runner.ran
    assert Potential.from_file(tmp_path / "task.pot").SCF_INFO.SCFSTATUS() == "CONVERGED"
