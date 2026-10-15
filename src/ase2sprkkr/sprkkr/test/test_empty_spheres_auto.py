from types import SimpleNamespace
from unittest.mock import patch

import pytest
from ase.build import bulk

from ..calculation import _resolve_empty_spheres
from ..calculator import SPRKKR


@pytest.mark.parametrize(
    ("input_parameters", "task", "should_add"),
    [
        ("SCF", None, True),
        ("DOS", None, False),
        ("SCF", "DOS", False),
        ("DOS", "SCF", True),
    ],
)
def test_auto_empty_spheres_only_for_effective_scf(
    tmp_path, input_parameters, task, should_add
):
    atoms = bulk("Fe")
    calculator = SPRKKR(atoms=atoms, directory=tmp_path, empty_spheres="auto")

    with patch("ase2sprkkr.sprkkr.calculation.add_empty_spheres") as add_empty_spheres:
        calculator.save_input(input_parameters=input_parameters, task=task)

    assert bool(add_empty_spheres.call_count) is should_add


def test_auto_empty_spheres_uses_canonical_execution_task_identity():
    parameters = SimpleNamespace(execution_task_name="SCF")
    atoms = SimpleNamespace(has_potential=lambda: False, sites=[])

    with patch("ase2sprkkr.sprkkr.calculation.add_empty_spheres") as add_empty_spheres:
        _resolve_empty_spheres(
            input_parameters=parameters,
            atoms=atoms,
            empty_spheres="auto",
            default_empty_spheres=False,
        )

    add_empty_spheres.assert_called_once_with(atoms)
