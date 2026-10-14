from ase.build import bulk
import pytest

from .. import es_finder


def test_radii_ratios_map_is_copied():
    original = dict(es_finder.winger_seitz_radii)

    radii = es_finder._prepare_radii_ratios_map(None)
    radii["DefinitelyNotAnElement"] = 1.0

    assert es_finder.winger_seitz_radii == original
    assert radii is not es_finder.winger_seitz_radii


def test_radii_ratios_map_overrides_defaults():
    radii = es_finder._prepare_radii_ratios_map({"Co": 123.0})

    assert radii["Co"] == 123.0
    assert radii.keys() == es_finder.winger_seitz_radii.keys()


@pytest.mark.skipif(not es_finder.is_enabled, reason="ES finder cannot be imported")
def test_es_finder():
    atoms = bulk("Co")
    assert len(es_finder.empty_spheres(atoms).positions) == 0

    atoms = bulk("NaCl", "rocksalt", a=5.64)
    out = es_finder.empty_spheres(atoms)
    assert len(out.positions) == 12
