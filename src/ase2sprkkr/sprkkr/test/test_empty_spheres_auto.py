from unittest.mock import patch

import pytest
from ase.build import bulk

from ..calculator import SPRKKR


@pytest.mark.parametrize(
    ("input_parameters", "should_add"),
    [
        ("SCF", True),
        ("DOS", False),
    ],
)
def test_auto_empty_spheres_only_for_scf(tmp_path, input_parameters, should_add):
    atoms = bulk("Fe")
    calculator = SPRKKR(atoms=atoms, directory=tmp_path, empty_spheres="auto")

    with patch("ase2sprkkr.sprkkr.calculation.add_empty_spheres") as add_empty_spheres:
        calculator.save_input(input_parameters=input_parameters)

    assert bool(add_empty_spheres.call_count) is should_add
