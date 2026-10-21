from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import Mock

import pytest
from ase import Atoms

from ase2sprkkr.common.file_fingerprints import file_fingerprint
from ase2sprkkr.potentials import potentials as potentials_module
from ase2sprkkr.potentials.potentials import Potential, PotentialFileChanged
from ase2sprkkr.sprkkr.configuration import ConfigurationFile


def _example_potential():
    return Path(__file__).parent.parent / "examples" / "FePt.new.pot"


def test_potential_attachment_sharing_and_reassignment():
    attached = Potential.from_file(_example_potential())
    atoms = attached.atoms

    assert atoms is not None
    assert atoms.has_potential()
    assert atoms.potential is attached

    # A second potential may use these atoms without replacing their owner.
    second = Potential.from_file(_example_potential(), atoms=atoms)
    assert second.atoms is atoms
    assert attached.atoms is atoms
    assert atoms.potential is attached

    # Reassigning a separate potential must not detach the first one.
    new_atoms = atoms.copy()
    standalone = Potential(atoms=atoms, set_from_atoms=False)
    standalone.atoms = new_atoms

    assert standalone.atoms is new_atoms
    assert attached.atoms is atoms
    assert atoms.potential is attached
    assert not new_atoms.has_potential()


def test_empty_atoms_are_explicit_set_from_atoms_input(monkeypatch):
    attached = Atoms("Fe")
    empty = Atoms()
    potential = Potential(atoms=attached, set_from_atoms=False)
    promoted = []

    monkeypatch.setattr(
        potentials_module.SPRKKRAtoms,
        "promote_ase_atoms",
        lambda atoms: promoted.append(atoms) or atoms,
    )
    monkeypatch.setattr(Atoms, "update_sites_from_symbols", lambda self: None, raising=False)
    monkeypatch.setattr(ConfigurationFile, "set_from_atoms", lambda self, atoms, io_data: None)

    potential.set_from_atoms(empty)

    assert promoted == [empty]


def test_save_to_file_synchronizes_atoms_once(monkeypatch):
    potential = Potential(atoms=Atoms("Fe"), set_from_atoms=False)
    empty = Atoms()
    calls = []

    monkeypatch.setattr(
        Potential,
        "set_from_atoms",
        lambda self, atoms=None, io_data=None: calls.append(atoms),
    )
    monkeypatch.setattr(
        ConfigurationFile,
        "save_to_file",
        lambda self, file, atoms=None, *, validate="save": file,
    )

    output = io.StringIO()
    potential.save_to_file(output)
    assert calls == [None]

    calls.clear()
    potential.save_to_file(output, atoms=empty)
    assert calls == [empty]

    calls.clear()
    potential.save_to_file(output, atoms=False)
    assert calls == [None]


def test_file_unchanged_check_is_opt_in(monkeypatch, tmp_path):
    # No real potential parsing is needed to test the provenance policy.
    path = tmp_path / "not-parsed.pot"
    parsed = object()
    reader = Mock(return_value=parsed)
    monkeypatch.setattr(Potential.potential_definition, "read_from_file", reader)
    monkeypatch.setattr(
        potentials_module,
        "file_fingerprint",
        lambda _path: pytest.fail("fingerprinting must be disabled by default"),
    )

    assert Potential.from_file(path) is parsed
    reader.assert_called_once_with(path, atoms=None, allow_dangerous=False)


def test_file_unchanged_check_rejects_rewrite(monkeypatch, tmp_path):
    # Fingerprint mismatch is checked around the parser, not by the parser.
    path = tmp_path / "not-parsed.pot"
    path.write_bytes(b"fake potential")
    original = file_fingerprint(path)
    changed = original._replace(size=original.size + 1)
    fingerprints = iter((original, changed))
    reader = Mock(return_value=object())
    monkeypatch.setattr(Potential.potential_definition, "read_from_file", reader)
    monkeypatch.setattr(
        potentials_module,
        "file_fingerprint",
        lambda _path: next(fingerprints),
    )

    with pytest.raises(PotentialFileChanged, match="changed while it was being read"):
        Potential.from_file(path, check_file_unchanged=True)
    reader.assert_called_once_with(path, atoms=None, allow_dangerous=False)
