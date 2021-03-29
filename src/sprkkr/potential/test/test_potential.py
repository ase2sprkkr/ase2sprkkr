if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..potentials import Potential
from ase.spacegroup import crystal
import io

class TestPotential(TestCase):

  def test_potential(self):
    a = 5.64
    nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
                   cellpar=[a, a, a, 90, 90, 90])

    sio = io.StringIO()
    p1 = Potential.from_atoms(nacl)
    p1.atoms.symbols[0] = 'K'
    p1.save_to_file(sio)
    sio.seek(0)
    p2 = Potential.from_file(sio)

    self.assertIsNotNone(p2.atoms)
    a1 = p1.atoms
    a2 = p2.atoms
    self.assertEqual(a2.cell * 1, a1.cell * 1)
    self.assertEqual(a2.positions, a1.positions)
    self.assertEqual(str(a2.symbols), str(a1.symbols))
