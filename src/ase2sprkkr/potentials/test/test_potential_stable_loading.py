from io import StringIO
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


def test_file_unchanged_check_is_opt_in(monkeypatch):
    monkeypatch.setattr(
        potentials_module,
        "file_fingerprint",
        lambda _path: pytest.fail("fingerprinting must be disabled by default"),
    )

    Potential.from_file(_example_potential())


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


def test_file_unchanged_check_requires_path():
    with pytest.raises(TypeError, match="requires a filesystem path"):
        Potential.from_file(StringIO(""), check_file_unchanged=True)
