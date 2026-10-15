from ..sprkkr_atoms import SPRKKRAtoms
from ..sites import SiteType


def _single_site():
    atoms = SPRKKRAtoms("Na", symmetry=False)
    site = atoms.sites[0]
    return atoms, site


def _assert_canonical_occupation(atoms, index=0):
    site = atoms.sites[index]
    kinds = atoms.arrays["spacegroup_kinds"]
    kind = int(kinds[index])
    expected = site.occupation.to_dict(True)

    assert atoms.numbers[index] == site.occupation.primary_atomic_number
    assert atoms.info["occupancy"][str(kind)] == expected

    reconstructed = SiteType.creator(atoms)(index)
    assert reconstructed.occupation.to_dict(True) == expected
    assert reconstructed.occupation.primary_atomic_number == atoms.numbers[index]


def test_delitem_keeps_canonical_ase_state_in_sync():
    atoms, site = _single_site()
    site.occupation.set({"Na": 0.6, "Fe": 0.4})

    del site.occupation["Fe"]

    _assert_canonical_occupation(atoms)
    assert site.occupation.to_dict(True) == {"Na": 0.6, "X": 0.4}


def test_clean_keeps_canonical_ase_state_in_sync():
    atoms, site = _single_site()
    site.occupation.set({"Na": 1.0, "Fe": 0.0})

    site.occupation.clean()

    _assert_canonical_occupation(atoms)
    assert site.occupation.to_dict(True) == {"Na": 1.0}


def test_replace_type_keeps_canonical_ase_state_in_sync():
    atoms, site = _single_site()

    site.occupation.replace_type("Na", "Fe")

    _assert_canonical_occupation(atoms)
    assert site.occupation.to_dict(True) == {"Fe": 1.0}
    assert atoms.symbols[0] == "Fe"
