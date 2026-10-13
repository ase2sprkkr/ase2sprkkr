if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import numpy as np
from ase import Atoms

from ..atoms_region import AtomsRegion  # NOQA E402


class TestAtomsRegion(TestCase):
    def test_slice_setter_invalidates_cached_membership(self):
        atoms = Atoms("HHeLi")
        region = AtomsRegion("test", slice(0, 2), np.eye(3))
        region.set_atoms(atoms, add=False)

        np.testing.assert_array_equal(region.ids, [0, 1])
        self.assertEqual(region.set_of_ids, {0, 1})

        region.slice = slice(1, None)

        np.testing.assert_array_equal(region.ids, [1, 2])
        self.assertEqual(region.set_of_ids, {1, 2})

    def test_create_atoms_returns_region_atoms(self):
        atoms = Atoms(
            "HHeLi",
            cell=np.diag([2.0, 3.0, 4.0]),
            pbc=[True, False, True],
        )
        region = AtomsRegion(
            "test",
            slice(1, 3),
            np.diag([5.0, 6.0, 7.0]),
            [False, True, False],
        )
        region.set_atoms(atoms, add=False)

        created = region.create_atoms()

        self.assertIsInstance(created, Atoms)
        self.assertEqual(created.get_chemical_symbols(), ["He", "Li"])
        np.testing.assert_allclose(created.cell, region.cell)
        np.testing.assert_array_equal(created.pbc, region.pbc)

    def test_spacegroup_info_uses_region_api(self):
        expected = object()

        class SpacegroupInfo:
            def info_for_region(inner_self, region):
                inner_self.region = region
                return expected

        class Parent:
            spacegroup_info = SpacegroupInfo()

        region = AtomsRegion("test", slice(None), np.eye(3))
        region.atoms = Parent()

        self.assertIs(region.spacegroup_info, expected)
        self.assertIs(region.atoms.spacegroup_info.region, region)

    def test_copy_for_atoms_preserves_incomplete_cell_and_pbc(self):
        parent = Atoms(
            cell=np.diag([2.0, 3.0, 4.0]),
            pbc=[True, False, True],
        )
        region = AtomsRegion(
            "test",
            slice(0, 0),
            np.array(
                [
                    [7.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 9.0],
                ]
            ),
            [False, None, True],
            atoms=parent,
        )
        new_parent = Atoms(
            cell=np.diag([11.0, 12.0, 13.0]),
            pbc=[True, True, False],
        )

        copied = region.copy_for_atoms(new_parent)

        np.testing.assert_allclose(copied.incomplete_cell, region.incomplete_cell)
        self.assertEqual(list(copied.incomplete_pbc), [False, None, True])
        np.testing.assert_allclose(copied.cell[0], [7.0, 0.0, 0.0])
        np.testing.assert_allclose(copied.cell[1], new_parent.cell[1])
        np.testing.assert_allclose(copied.cell[2], [0.0, 0.0, 9.0])
        np.testing.assert_array_equal(copied.pbc, [False, True, True])
