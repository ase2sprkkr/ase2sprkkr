import os

import pytest
from ase import io
from ase.build import bulk

from ..occupations import Occupation
from ..sprkkr_atoms import SPRKKRAtoms


def test_occupation_updates_symbols_and_respects_symmetry():
    atoms = io.read(os.path.join(os.path.dirname(__file__), "test.cif"))
    SPRKKRAtoms.promote_ase_atoms(atoms)

    assert atoms.sites[0].occupation["Hf"] == 0.7025
    assert atoms.sites[0].occupation["Ti"] == 0.2975
    assert atoms.sites[-1].occupation["H"] == 0.233

    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    SPRKKRAtoms.promote_ase_atoms(atoms)
    with pytest.raises(RuntimeError):
        atoms.sites[0].occupation["Fe"] = 0.5

    atoms.sites[0].site_type.occupation["Fe"] = 0.4
    atoms.sites[0].site_type.occupation["Ni"] = 0.6
    assert atoms.symbols[0] == "Ni"

    atoms.sites[0].break_symmetry()
    atoms.sites[0].occupation["Fe"] = 0.7
    assert atoms.symbols[0] == "Fe"


def test_occupation_public_methods():
    occupation = Occupation({"Fe": 0.4, "Ni": 0.6})

    occupation.check()
    values = {atomic_type.symbol: value for atomic_type, value in occupation.to_tuple()}
    assert values == {"Fe": 0.4, "Ni": 0.6}

    occupation._occupation[occupation.atomic_type("Fe")] = 0.3
    with pytest.raises(ValueError, match="equal to one"):
        occupation.check()
