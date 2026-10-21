import numpy as np
import pytest
from ase.build import bulk

from ..build import semiinfinite_system, stack
from ..sprkkr_atoms import SPRKKRAtoms


@pytest.fixture(scope="module")
def cu_sc():
    return bulk("Cu", "sc", a=2.0)


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

    @pytest.mark.parametrize("orientation_arg", ("hkl", "hkl2"))
    def test_deprecated_orientation_builds_and_warns(
        self, monkeypatch, cu_sc, orientation_arg,
    ):
        # Keep the original end-to-end regression: the deprecated arguments
        # must still build a system, not just dispatch to rotate().
        from .. import build as build_module

        real_rotate = build_module.rotate
        rotations = []

        def record_rotate(atoms, direction):
            rotations.append((atoms, direction))
            return real_rotate(atoms, direction)

        monkeypatch.setattr(build_module, "rotate", record_rotate)
        atoms = cu_sc.copy()
        options = {orientation_arg: (0, 0, 1)}
        if orientation_arg == "hkl2":
            options["atoms2"] = atoms.copy()

        with pytest.warns(
            DeprecationWarning, match=rf"{orientation_arg} parameter",
        ):
            system = semiinfinite_system(atoms, 1, **options)

        assert len(system) > 0
        assert len(rotations) == 1
        assert rotations[0][0] is options.get("atoms2", atoms)
        assert rotations[0][1] == (0, 0, 1)

    def test_zero_repeat(self, cu_sc):
        atoms = cu_sc.copy()
        SPRKKRAtoms.promote_ase_atoms(atoms)

        # All three variants use the same pristine bulk structure. Construct
        # it and initialize its sites once, then exercise each repeat policy.
        for repeat, symbols in (
            ((0, 0), "CuX"),
            ((2, 0), "Cu3X"),
            ((0, 2), "CuX3"),
        ):
            system = semiinfinite_system(atoms, repeat)

            assert str(system.symbols) == symbols
            np.testing.assert_allclose(
                system.positions, (atoms * (1, 1, len(system))).positions
            )
            np.testing.assert_allclose(system.cell[2], len(system) * atoms.cell[2])
            np.testing.assert_array_equal(system.regions["left"].ids, [0])
            np.testing.assert_array_equal(
                system.regions["central"].ids, np.arange(1, len(system) - 1)
            )
            np.testing.assert_array_equal(
                system.regions["right"].ids, [len(system) - 1]
            )
            assert len(np.unique(system.positions, axis=0)) == len(system)

    def test_stack_last_region_is_open_ended(self, cu_sc):
        atoms = cu_sc.copy()
        system = stack({"first": atoms, "last": atoms}, axis=2)

        assert system.regions["first"].slice.stop == len(atoms)
        assert system.regions["last"].slice.stop is None
