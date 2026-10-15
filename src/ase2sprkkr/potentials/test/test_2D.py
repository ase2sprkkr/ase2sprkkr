from pathlib import Path

import numpy as np

from ..potentials import Potential


EXAMPLE_POTENTIAL = Path(__file__).parent.parent / "examples" / "GeTe.pot"


def test_2d_potential_round_trip():
    potential = Potential.from_file(EXAMPLE_POTENTIAL)
    atoms = potential.atoms

    assert str(atoms.symbols) == "GeXTeX" * 8 + "GeXTeX13"
    assert len(atoms.regions["left"]) == 4
    assert len(atoms.regions["right"]) == 4

    restored = Potential.from_string(potential.to_string()).atoms

    np.testing.assert_allclose(atoms.cell, restored.cell)
    np.testing.assert_allclose(atoms.regions["left"].cell, restored.regions["left"].cell)
    np.testing.assert_allclose(atoms.regions["right"].cell, restored.regions["right"].cell)
    np.testing.assert_allclose(atoms.positions, restored.positions)
    assert len(atoms.regions["left"]) == len(restored.regions["left"])
    assert len(atoms.regions["right"]) == len(restored.regions["right"])
    assert str(atoms.symbols) == str(restored.symbols)
