import numpy as np
from ase import Atoms
from ase.build import bulk

from ..sprkkr_atoms import SPRKKRAtoms


class TestSPRKKRAtoms:
    def test_extend(self):
        atoms = SPRKKRAtoms("NaCl")
        atoms.set_cell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        atoms.set_positions([[0, 0, 0], [0, 1, 0]])

        added = SPRKKRAtoms("Na")
        added.set_cell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        added.set_positions([[1, 0, 0]])

        original_site_type_id = id(atoms.sites[0].site_type)
        added_site_type_id = id(added.sites[0].site_type)
        atoms.extend(added)

        assert id(atoms.sites[0].site_type) == original_site_type_id
        assert id(atoms.sites[2].site_type) != added_site_type_id
        assert added.sites[0] is not atoms.sites[2]

        ase_atoms = Atoms("NaCl")
        ase_atoms.set_cell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        ase_atoms.set_positions([[0, 0, 0], [0, 1, 0]])
        ase_atoms.extend(added)
        SPRKKRAtoms.promote_ase_atoms(ase_atoms)
        assert id(ase_atoms.sites[2].site_type) != added_site_type_id

    def test_extend_spacegroup_kinds(self):
        atoms = SPRKKRAtoms("NaCl")
        atoms.arrays["spacegroup_kinds"] = np.array([0, 1])
        atoms.info["occupancy"] = {
            "0": {"Na": 1.0},
            "1": {"Cl": 1.0},
        }

        added = SPRKKRAtoms("NaCl")
        added.arrays["spacegroup_kinds"] = np.array([0, 1])
        added.info["occupancy"] = {
            "0": {"Na": 0.5, "K": 0.5},
            "1": {"Cl": 0.5, "I": 0.5},
        }

        atoms.extend(added)

        np.testing.assert_array_equal(atoms.arrays["spacegroup_kinds"], [0, 1, 2, 3])
        assert atoms.info["occupancy"]["0"] == {"Na": 1.0}
        assert atoms.info["occupancy"]["1"] == {"Cl": 1.0}
        assert atoms.info["occupancy"]["2"] == {"Na": 0.5, "K": 0.5}
        assert atoms.info["occupancy"]["3"] == {"Cl": 0.5, "I": 0.5}

        atoms = SPRKKRAtoms("Na")
        atoms.arrays["spacegroup_kinds"] = np.array([0])
        added = SPRKKRAtoms("Cl")
        added.arrays["spacegroup_kinds"] = np.array([0])

        atoms.extend(added)

        np.testing.assert_array_equal(atoms.arrays["spacegroup_kinds"], [0, 1])

    def test_atoms(self):
        atoms = Atoms("NaCl")
        atoms.set_cell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        atoms.set_positions([[0, 0, 0], [0, 1, 0]])
        atoms.info["occupancy"] = {0: {"Na": 1.0}, 1: {"Cl": 0.4, "I": 0.6}}
        atoms.arrays["spacegroup_kinds"] = np.asarray([0, 1])

        assert str(atoms.symbols) == "NaCl"
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert len(atoms.sites) == 2
        assert len(atoms.sites[0].occupation) == 1
        assert len(atoms.sites[1].occupation) == 2
        assert atoms.sites[1].occupation["Cl"] == 0.4
        assert str(atoms.symbols) == "NaI"
        assert atoms.sites[0].occupation.as_dict == {"Na": 1}
        assert atoms.sites[1].occupation.as_dict == {"Cl": 0.4, "I": 0.6}
        assert atoms.info["occupancy"]["1"] == {"Cl": 0.4, "I": 0.6}

        atoms.sites[1].site_type = atoms.sites[0].site_type
        atoms.compute_sites_symmetry()

    def test_symmetry(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert atoms.sites[1].site_type == atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.symmetry = False
        assert atoms.sites[1].site_type != atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.symmetry = True
        assert atoms.sites[1].site_type == atoms.sites[3].site_type
        atoms.sites[3].site_type = atoms.sites[3].site_type.copy()
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.symmetry = True
        assert atoms.sites[1].site_type != atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.compute_sites_symmetry()
        assert atoms.sites[1].site_type == atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.symmetry = False
        assert atoms.sites[1].site_type != atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.sites[3].site_type = atoms.sites[1].site_type
        atoms.symmetry = False
        assert atoms.sites[1].site_type == atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms.breaks_sites_symmetry()
        assert atoms.sites[1].site_type != atoms.sites[3].site_type
        assert atoms.sites[1] is not atoms.sites[3]

        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        SPRKKRAtoms.promote_ase_atoms(atoms, symmetry=False)
        assert atoms.sites[1].site_type != atoms.sites[3].site_type

        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert atoms.sites[1].site_type != atoms.sites[3].site_type

        SPRKKRAtoms.promote_ase_atoms(atoms, symmetry=True)
        assert atoms.sites[1].site_type == atoms.sites[3].site_type

    def test_occupancy(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert atoms.sites[1].site_type == atoms.sites[3].site_type
        atoms.breaks_sites_symmetry()
        assert atoms.sites[1].site_type != atoms.sites[3].site_type

        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        atoms.info["occupancy"] = {
            0: {"Li": 1},
            1: {"Cl": 1},
            2: {"Li": 1},
            3: {"Cl": 1},
        }
        atoms.arrays["spacegroup_kinds"] = np.array([0, 1, 0, 1])
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert atoms.sites[1].site_type == atoms.sites[3].site_type

        atoms = bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
        atoms.info["occupancy"] = {
            0: {"Li": 1},
            1: {"Cl": 0.5, "I": 0.5},
            2: {"Li": 1},
            3: {"Cl": 1},
        }
        atoms.arrays["spacegroup_kinds"] = np.array([0, 1, 0, 2])
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert atoms.sites[1].site_type != atoms.sites[3].site_type
        assert atoms.sites[1].occupation.as_dict == {"Cl": 0.5, "I": 0.5}
