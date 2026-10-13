from ase.build import bulk
import pytest

from .. import es_finder


@pytest.mark.skipif(not es_finder.is_enabled, reason="ES finder cannot be imported")
def test_es_finder():
    atoms = bulk("Co")
    assert len(es_finder.empty_spheres(atoms).positions) == 0

    atoms = bulk("NaCl", "rocksalt", a=5.64)
    out = es_finder.empty_spheres(atoms)
    assert len(out.positions) == 12
