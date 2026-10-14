from __future__ import annotations

from pathlib import Path

import pytest

from ase2sprkkr.common.file_fingerprints import file_fingerprint
from ase2sprkkr.potentials import potentials as potentials_module
from ase2sprkkr.potentials.potentials import Potential, PotentialFileChanged


def _example_potential():
    return Path(__file__).parent.parent / "examples" / "FePt.new.pot"


def test_potential_from_file_attaches_itself_to_loaded_atoms():
    potential = Potential.from_file(_example_potential())
    atoms = potential.atoms

    assert atoms is not None
    assert atoms.has_potential()
    assert atoms.potential is potential


def test_potential_from_file_attaches_to_supplied_atoms():
    source = Potential.from_file(_example_potential())
    atoms = source.atoms.copy()
    assert not atoms.has_potential()

    potential = Potential.from_file(_example_potential(), atoms=atoms)

    assert potential.atoms is atoms
    assert atoms.potential is potential


def test_replacing_attached_potential_detaches_previous_one():
    previous = Potential.from_file(_example_potential())
    atoms = previous.atoms

    replacement = Potential.from_file(_example_potential(), atoms=atoms)

    assert atoms.potential is replacement
    assert replacement.atoms is atoms
    assert previous._atoms is None

    previous_atoms = previous.atoms
    assert previous_atoms is not atoms
    assert previous_atoms.potential is previous


def test_reassigning_potential_atoms_detaches_old_atoms():
    potential = Potential.from_file(_example_potential())
    old_atoms = potential.atoms
    new_atoms = old_atoms.copy()

    potential.atoms = new_atoms

    assert not old_atoms.has_potential()
    assert new_atoms.potential is potential
    assert potential.atoms is new_atoms


def test_file_unchanged_check_is_opt_in(monkeypatch):
    path = _example_potential()

    monkeypatch.setattr(
        potentials_module,
        "file_fingerprint",
        lambda _path: pytest.fail("fingerprinting must be disabled by default"),
    )

    Potential.from_file(path)


def test_file_unchanged_check_rejects_rewrite(monkeypatch):
    path = _example_potential()
    original = file_fingerprint(path)
    changed = original._replace(size=original.size + 1)
    fingerprints = iter((original, changed))

    monkeypatch.setattr(
        potentials_module,
        "file_fingerprint",
        lambda _path: next(fingerprints),
    )

    with pytest.raises(PotentialFileChanged, match="changed while it was being read"):
        Potential.from_file(path, check_file_unchanged=True)
