from ase.build import bulk
import numpy as np
import pytest

from ..build import aperiodic_times, stack, rotate, _minimal_surface_layers


@pytest.fixture(scope="module")
def licl_rocksalt():
    return bulk("LiCl", "rocksalt", a=5.64)


@pytest.fixture(scope="module")
def cu_sc():
    return bulk("Cu", "sc", a=2.0)


class TestBuild:
    def test(self, licl_rocksalt):
        atoms = licl_rocksalt.copy()
        assert len(aperiodic_times(atoms, (2.1, 1, 1))) == 6
        assert len(aperiodic_times(atoms, (1, 2.1, 1))) == 5
        assert len(aperiodic_times(atoms, (1, 2.9, 1), direction=-1)) == 5
        at = aperiodic_times(atoms, (1, 3.1, 1))
        assert len(at) == 7
        assert len(aperiodic_times(atoms, (1, 3.1, 2.1))) == 18
        np.testing.assert_allclose(at.positions[::2], np.arange(4)[:, None] * atoms.cell[1] + atoms.positions[0])
        np.testing.assert_allclose(at.positions[1::2], np.arange(3)[:, None] * atoms.cell[1] + atoms.positions[1])
        assert str(at.symbols) == "LiClLiClLiClLi"
        np.testing.assert_array_equal(at.pbc, [True, False, True])
        cell = atoms.cell.copy()
        cell[1] *= 3.1
        np.testing.assert_allclose(at.cell, cell)

        at = aperiodic_times(atoms, (1, 3.9, 1), direction=-1)
        np.testing.assert_allclose(
            at.positions[::2], np.arange(4)[:, None] * atoms.cell[1] + atoms.positions[1] - 0.1 * atoms.cell[1]
        )
        np.testing.assert_allclose(
            at.positions[1::2], np.arange(1, 4)[:, None] * atoms.cell[1] + atoms.positions[0] - 0.1 * atoms.cell[1]
        )
        assert str(at.symbols) == "ClLiClLiClLiCl"
        np.testing.assert_array_equal(at.pbc, [True, False, True])
        cell[1] = atoms.cell[1] * 3.9
        np.testing.assert_allclose(at.cell, cell)

    def test_stack(self, licl_rocksalt):
        atoms = licl_rocksalt.copy()

        np.testing.assert_allclose(stack([atoms, atoms, atoms], axis=2).positions, (atoms * (1, 1, 3)).positions)

        a2 = atoms.copy()
        a2.symbols = "KN"
        a2.positions = a2.positions + 1
        stacked = stack([atoms, a2], axis=2)
        assert str(stacked.symbols) == "LiClKN"
        np.testing.assert_allclose(stacked.positions, np.concatenate([atoms.positions, a2.positions + atoms.cell[2]]))

    def test_stack_scalar_at_uses_stacking_cell_vector_length(self, cu_sc):
        atoms = cu_sc.copy()

        stacked = stack([atoms, atoms], axis=2, at=[None, 3.0, 6.0])

        np.testing.assert_allclose(stacked.positions[:, 2], [0.0, 3.0])
        np.testing.assert_allclose(stacked.cell[2], [0.0, 0.0, 6.0])

    def test_stack_accepts_iterator(self, cu_sc):
        atoms = cu_sc.copy()

        stacked = stack((atoms.copy() for _ in range(3)), axis=2)
        expected = atoms * (1, 1, 3)

        np.testing.assert_allclose(stacked.positions, expected.positions)
        np.testing.assert_allclose(stacked.cell, expected.cell)

    def test_stack_periodic_argument_wins_on_stacking_axis(self, cu_sc):
        atoms = cu_sc.copy()
        atoms.pbc[2] = False

        stacked = stack([atoms, atoms], axis=2, periodic=True)

        assert stacked.pbc[2]

    def test_stack_empty_atoms(self, cu_sc):
        atoms = cu_sc.copy()
        empty = atoms * (1, 1, 0)
        expected = atoms * (1, 1, 2)

        middle_empty = stack([atoms, empty, atoms], axis=2)
        np.testing.assert_allclose(middle_empty.positions, expected.positions)
        np.testing.assert_allclose(middle_empty.cell, expected.cell)

        trailing_empty = stack([atoms, atoms, empty], axis=2)
        np.testing.assert_allclose(trailing_empty.positions, expected.positions)
        np.testing.assert_allclose(trailing_empty.cell[2], 2 * atoms.cell[2])

    def test_minimal_surface_layers(self):
        atoms = bulk("Cu", "fcc", cubic=True)

        assert _minimal_surface_layers(atoms, (1, 0, 0)) == 1
        assert _minimal_surface_layers(atoms, (1, 1, 0)) == 2
        assert _minimal_surface_layers(atoms, (1, 1, 1)) == 3

    def test_rotate_fcc_au_111(self):
        atoms = bulk("Au", "fcc", cubic=True)
        rotated = rotate(atoms, (1, 1, 1))

        assert _minimal_surface_layers(atoms, (1, 1, 1)) == 3
        np.testing.assert_array_equal(rotated.pbc, [True, True, True])
        assert len(rotated) == 3 * len(atoms)
