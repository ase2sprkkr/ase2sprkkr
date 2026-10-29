from types import SimpleNamespace
from unittest.mock import patch

import pytest
from ase.build import bulk

from .. import _preparation as preparation
from .._preparation import resolve_empty_spheres
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
    tmp_path, monkeypatch, input_parameters, task, should_add
):
    atoms = bulk("Fe")
    calculator = SPRKKR(atoms=atoms, directory=tmp_path, empty_spheres="auto")

    class PolicyChecked(Exception):
        """Abort after checking the real policy, before file serialization."""

    selected_tasks = []

    def check_policy(**kwargs):
        selected_tasks.append(kwargs["input_parameters"].task_name)
        resolve_empty_spheres(**kwargs)
        raise PolicyChecked

    monkeypatch.setattr(preparation, "resolve_empty_spheres", check_policy)
    with patch("ase2sprkkr.sprkkr._preparation.add_empty_spheres") as add_empty_spheres:
        with pytest.raises(PolicyChecked):
            calculator.save_input(input_parameters=input_parameters, task=task)

    assert selected_tasks == [task or input_parameters]
    assert bool(add_empty_spheres.call_count) is should_add
    assert not list(tmp_path.glob("*.inp"))
    assert not list(tmp_path.glob("*.pot"))


def test_auto_empty_spheres_uses_effective_task_identity():
    parameters = SimpleNamespace(task_name="SCF")
    atoms = SimpleNamespace(has_potential=lambda: False, sites=[])

    with patch("ase2sprkkr.sprkkr._preparation.add_empty_spheres") as add_empty_spheres:
        resolve_empty_spheres(
            input_parameters=parameters,
            atoms=atoms,
            empty_spheres="auto",
            default_empty_spheres=False,
        )

    add_empty_spheres.assert_called_once_with(atoms)
