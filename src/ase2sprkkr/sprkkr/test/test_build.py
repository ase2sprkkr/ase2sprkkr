if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import pytest
import numpy as np

from ..build import semiinfinite_system, stack  # NOQA E402
from ..sprkkr_atoms import SPRKKRAtoms  # NOQA E402
from ase.build import bulk  # NOQA E402


class TestBuild(TestCase):
    @pytest.mark.slow
    def test(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64)
        si = semiinfinite_system(atoms, (4.2))
        self.assertEqual(len(si), 22)
        self.assertEqual(sum(si.symbols == "Li"), 6)
        self.assertEqual(sum(si.symbols == "Cl"), 5)
        self.assertEqual(sum(si.symbols == "X"), 11)
        self.assertEqual(si.positions, (atoms * (1, 1, 11)).positions)

    def test_deprecated_hkl_warns(self):
        atoms = bulk("Cu", "sc", a=2.0)
        with pytest.warns(DeprecationWarning, match="hkl parameter"):
            semiinfinite_system(atoms, 1, hkl=(0, 0, 1))

    def test_deprecated_hkl2_warns(self):
        atoms = bulk("Cu", "sc", a=2.0)
        with pytest.warns(DeprecationWarning, match="hkl2 parameter"):
            semiinfinite_system(atoms, 1, atoms2=atoms.copy(), hkl2=(0, 0, 1))

    @pytest.mark.parametrize(
        "repeat, symbols",
        (
            ((0, 0), "CuX"),
            ((2, 0), "Cu3X"),
            ((0, 2), "CuX3"),
        ),
    )
    def test_zero_repeat(self, repeat, symbols):
        atoms = bulk("Cu", "sc", a=2.0)
        SPRKKRAtoms.promote_ase_atoms(atoms)
        system = semiinfinite_system(atoms, repeat)

        self.assertEqual(str(system.symbols), symbols)
        self.assertEqual(system.positions, (atoms * (1, 1, len(system))).positions)
        self.assertEqual(system.cell[2], len(system) * atoms.cell[2])
        self.assertEqual(system.regions["left"].ids, np.array([0]))
        self.assertEqual(system.regions["central"].ids, np.arange(1, len(system) - 1))
        self.assertEqual(system.regions["right"].ids, np.array([len(system) - 1]))
        self.assertEqual(len(np.unique(system.positions, axis=0)), len(system))

    def test_stack_last_region_is_open_ended(self):
        atoms = bulk("Cu", "sc", a=2.0)
        system = stack({"first": atoms, "last": atoms}, axis=2)

        self.assertEqual(system.regions["first"].slice.stop, len(atoms))
        self.assertIsNone(system.regions["last"].slice.stop)
