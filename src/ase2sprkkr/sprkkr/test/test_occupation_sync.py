import numpy as np
from ase.build import bulk

from ..sprkkr_atoms import SPRKKRAtoms


def _kind_occupancy(atoms, index=0):
    kind = int(atoms.arrays["spacegroup_kinds"][index])
    return atoms.info["occupancy"][str(kind)]


def test_occupation_mutations_keep_ase_state_synchronized():
    atoms = SPRKKRAtoms("Na")
    site = atoms.sites[0]

    site.occupation = {"Na": 0.6, "K": 0.4}
    assert atoms.symbols[0] == "Na"
    assert _kind_occupancy(atoms) == {"Na": 0.6, "K": 0.4}

    site.occupation.set({"Na": 0.7, "K": 0.3})
    assert atoms.symbols[0] == "Na"
    assert _kind_occupancy(atoms) == {"Na": 0.7, "K": 0.3}

    site.occupation.replace_type("Na", "Cl")
    assert atoms.symbols[0] == "Cl"
    assert _kind_occupancy(atoms) == {"Cl": 0.7, "K": 0.3}

    site.occupation["Cl"] = 0.0
    site.occupation.clean()
    assert atoms.symbols[0] == "K"
    assert _kind_occupancy(atoms) == {"K": 1.0}

    site.occupation.add("Na", 0.25)
    assert atoms.symbols[0] == "K"
    assert _kind_occupancy(atoms) == {"K": 0.75, "Na": 0.25}

    del site.occupation["K"]
    assert atoms.symbols[0] == "Na"
    assert _kind_occupancy(atoms) == {"Na": 0.25, "X": 0.75}


def test_site_type_occupation_update_is_incremental(monkeypatch):
    atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
    SPRKKRAtoms.promote_ase_atoms(atoms)
    sites = atoms.sites

    site_type = sites[1].site_type
    indexes = site_type.index()
    kinds_before = atoms.arrays["spacegroup_kinds"].copy()
    kind = int(kinds_before[indexes[0]])
    assert np.all(kinds_before[indexes] == kind)

    def full_rebuild_must_not_run(*args, **kwargs):
        raise AssertionError("ordinary occupation mutation rebuilt spacegroup kinds")

    monkeypatch.setattr(atoms.spacegroup_info, "update_spacegroup_kinds", full_rebuild_must_not_run)

    site_type.occupation.replace_type("Cl", "I")

    np.testing.assert_array_equal(atoms.arrays["spacegroup_kinds"], kinds_before)
    assert np.all(atoms.numbers[indexes] == 53)
    assert atoms.info["occupancy"][str(kind)] == {"I": 1.0}


def test_site_type_update_rebuilds_missing_kind_mapping(monkeypatch):
    atoms = SPRKKRAtoms("Na")
    site_type = atoms.sites[0].site_type
    del atoms.arrays["spacegroup_kinds"]

    called = []
    original = atoms.spacegroup_info.update_spacegroup_kinds

    def record_rebuild(*args, **kwargs):
        called.append(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(atoms.spacegroup_info, "update_spacegroup_kinds", record_rebuild)

    site_type.occupation.replace_type("Na", "Cl")

    assert called == [{"if_required": True}]
    assert atoms.symbols[0] == "Cl"
    assert _kind_occupancy(atoms) == {"Cl": 1.0}


def test_occupation_reader_uses_spacegroup_kind_not_atom_index():
    atoms = SPRKKRAtoms(
        "Na2",
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        cell=[3.0, 3.0, 3.0],
        pbc=True,
        symmetry=False,
    )
    atoms.arrays["spacegroup_kinds"] = np.array([7, 7])
    atoms.info["occupancy"] = {"7": {"Na": 0.25, "K": 0.75}}

    sites = atoms.sites

    assert sites[0].occupation.as_dict == {"Na": 0.25, "K": 0.75}
    assert sites[1].occupation.as_dict == {"Na": 0.25, "K": 0.75}
