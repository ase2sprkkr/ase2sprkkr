from ase.build import bulk
import numpy as np

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..build import aperiodic_times, stack, rotate, _minimal_surface_layers  # NOQA: E402


class TestBuild(TestCase):
    def test(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64)
        assert len(aperiodic_times(atoms, (2.1, 1, 1))) == 6
        assert len(aperiodic_times(atoms, (1, 2.1, 1))) == 5
        assert len(aperiodic_times(atoms, (1, 2.9, 1), direction=-1)) == 5
        assert len(aperiodic_times(atoms, (1, 3.1, 1))) == 7
        assert len(aperiodic_times(atoms, (1, 3.1, 2.1))) == 18
        at = aperiodic_times(atoms, (1, 3.1, 1))
        self.assertEqual(at.positions[::2], np.arange(4)[:, None] * atoms.cell[1] + atoms.positions[0])
        self.assertEqual(at.positions[1::2], np.arange(3)[:, None] * atoms.cell[1] + atoms.positions[1])
        self.assertEqual("LiClLiClLiClLi", str(at.symbols))
        self.assertEqual(at.pbc, np.asarray([True, False, True]))
        cell = atoms.cell.copy()
        cell[1] *= 3.1
        assert (at.cell == cell).all()

        at = aperiodic_times(atoms, (1, 3.9, 1), direction=-1)
        self.assertEqual(
            at.positions[::2], np.arange(4)[:, None] * atoms.cell[1] + atoms.positions[1] - 0.1 * atoms.cell[1]
        )
        self.assertEqual(
            at.positions[1::2], np.arange(1, 4)[:, None] * atoms.cell[1] + atoms.positions[0] - 0.1 * atoms.cell[1]
        )
        self.assertEqual("ClLiClLiClLiCl", str(at.symbols))
        self.assertEqual(at.pbc, np.asarray([True, False, True]))
        cell[1] = atoms.cell[1] * 3.9
        assert (at.cell == cell).all()

    def test_stack(self):
        atoms = bulk("LiCl", "rocksalt", a=5.64)

        self.assertEqual(stack([atoms, atoms, atoms], axis=2).positions, (atoms * (1, 1, 3)).positions)

        a2 = atoms.copy()
        a2.symbols = "KN"
        a2.positions = a2.positions + 1
        self.assertEqual(str(stack([atoms, a2], axis=2).symbols), "LiClKN")
        self.assertEqual(
            stack([atoms, a2], axis=2).positions, np.concatenate([atoms.positions, a2.positions + atoms.cell[2]])
        )

    def test_minimal_surface_layers(self):
        atoms = bulk("Cu", "fcc", cubic=True)

        self.assertEqual(_minimal_surface_layers(atoms, (1, 0, 0)), 1)
        self.assertEqual(_minimal_surface_layers(atoms, (1, 1, 0)), 2)
        self.assertEqual(_minimal_surface_layers(atoms, (1, 1, 1)), 3)

    def test_rotate_fcc_au_111(self):
        atoms = bulk("Au", "fcc", cubic=True)
        rotated = rotate(atoms, (1, 1, 1))

        self.assertEqual(_minimal_surface_layers(atoms, (1, 1, 1)), 3)
        self.assertEqual(rotated.pbc, np.asarray([True, True, True]))
        self.assertEqual(len(rotated), 3 * len(atoms))
        self.assertEqual(np.unique(rotated.get_tags()), np.array([1, 2, 3]))
