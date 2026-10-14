import numpy as np
import pytest
from ase.build import bulk

from ..build import semiinfinite_system, stack
from ..sprkkr_atoms import SPRKKRAtoms


class TestBuild:
    @pytest.mark.slow
    def test(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64)
        si = semiinfinite_system(atoms, 4.2)
        assert len(si) == 22
        assert sum(si.symbols == "Li") == 6
        assert sum(si.symbols == "Cl") == 5
        assert sum(si.symbols == "X") == 11
        np.testing.assert_allclose(si.positions, (atoms * (1, 1, 11)).positions)

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

        assert str(system.symbols) == symbols
        np.testing.assert_allclose(system.positions, (atoms * (1, 1, len(system))).positions)
        np.testing.assert_allclose(system.cell[2], len(system) * atoms.cell[2])
        np.testing.assert_array_equal(system.regions["left"].ids, [0])
        np.testing.assert_array_equal(system.regions["central"].ids, np.arange(1, len(system) - 1))
        np.testing.assert_array_equal(system.regions["right"].ids, [len(system) - 1])
        assert len(np.unique(system.positions, axis=0)) == len(system)

    def test_stack_last_region_is_open_ended(self):
        atoms = bulk("Cu", "sc", a=2.0)
        system = stack({"first": atoms, "last": atoms}, axis=2)

        assert system.regions["first"].slice.stop == len(atoms)
        assert system.regions["last"].slice.stop is None
