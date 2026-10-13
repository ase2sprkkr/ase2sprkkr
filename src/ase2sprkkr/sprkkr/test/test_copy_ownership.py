import numpy as np

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..sprkkr_atoms import SPRKKRAtoms  # NOQA: E402
from ..sites import Site  # NOQA: E402


class TestSPRKKRAtomsCopyOwnership(TestCase):
    def test_copy_sites_uses_target_atoms(self):
        atoms = SPRKKRAtoms("Fe2", positions=[[0, 0, 0], [1, 0, 0]], symmetry=False)
        atoms.sites[1].site_type = atoms.sites[0].site_type
        target = SPRKKRAtoms("Fe2", positions=atoms.positions, symmetry=False)

        copied_sites = Site.copy_sites(atoms.sites, target)

        self.assertIs(copied_sites[0].site_type, copied_sites[1].site_type)
        self.assertIs(copied_sites[0].site_type.atoms, target)
        self.assertIsNot(copied_sites[0].site_type, atoms.sites[0].site_type)

    def test_copy_owns_independent_site_state(self):
        atoms = SPRKKRAtoms("Fe", symmetry=False)
        site = atoms.sites[0]
        site_type = site.site_type
        atomic_type = site_type.atomic_type(0)

        atomic_type.potential = np.ones((2, len(site_type.mesh)))
        atomic_type.charge = np.ones(len(site_type.mesh))
        atomic_type.moments = (1.0, 2.0, 3.0, 4.0, 5.0)
        atomic_type.magnetisation[:] = [1.0, 2.0, 3.0]
        site.magnetisation[:] = [4.0, 5.0]

        copied = atoms.copy()
        copied_site = copied.sites[0]
        copied_site_type = copied_site.site_type
        copied_atomic_type = copied_site_type.atomic_type(0)

        self.assertIs(copied_site_type.atoms, copied)
        self.assertIsNot(copied_site_type, site_type)
        self.assertIsNot(copied_site_type.mesh, site_type.mesh)
        self.assertIsNot(copied_atomic_type, atomic_type)
        self.assertIs(copied_atomic_type.mesh, copied_site_type.mesh)
        self.assertIsNot(copied_atomic_type.potential, atomic_type.potential)
        self.assertIsNot(copied_atomic_type.charge, atomic_type.charge)
        self.assertIsNot(copied_atomic_type.moments, atomic_type.moments)
        self.assertIs(copied_atomic_type.potential.mesh, copied_site_type.mesh)
        self.assertIs(copied_atomic_type.charge.mesh, copied_site_type.mesh)
        self.assertIsNot(copied_atomic_type.magnetisation, atomic_type.magnetisation)
        self.assertIsNot(copied_site.magnetisation, site.magnetisation)

        copied_atomic_type.potential.raw_value[0, 0] = 10.0
        copied_atomic_type.charge.raw_value[0] = 20.0
        copied_atomic_type.moments.smt = 30.0
        copied_atomic_type.magnetisation[0] = 40.0
        copied_site.magnetisation[0] = 50.0

        self.assertEqual(atomic_type.potential.raw_value[0, 0], 1.0)
        self.assertEqual(atomic_type.charge.raw_value[0], 1.0)
        self.assertEqual(atomic_type.moments.smt, 3.0)
        self.assertEqual(atomic_type.magnetisation[0], 1.0)
        self.assertEqual(site.magnetisation[0], 4.0)

    def test_copy_preserves_symmetry_sharing_only_inside_copy(self):
        atoms = SPRKKRAtoms("Fe2", positions=[[0, 0, 0], [1, 0, 0]], symmetry=False)
        atoms.sites[1].site_type = atoms.sites[0].site_type

        copied = atoms.copy()

        self.assertIs(copied.sites[0].site_type, copied.sites[1].site_type)
        self.assertIsNot(copied.sites[0].site_type, atoms.sites[0].site_type)
        self.assertIs(copied.sites[0].site_type.atoms, copied)
