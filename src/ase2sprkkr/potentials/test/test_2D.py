import os

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..potentials import Potential   # NOQA E402


class Test2DPotential(TestCase):

    def test(self):
        p=Potential.from_file(os.path.join(os.path.dirname(__file__),'..','examples', 'GeTe.pot'))
        atoms = p.atoms
        self.assertEqual(str(atoms.symbols), "GeXTeX" * 8 + "GeXTeX13")
        self.assertEqual(len(atoms.regions['left']), 4)
        self.assertEqual(len(atoms.regions['right']), 4)

        pf = p.to_string()
        p2=Potential.from_string(pf)
        a2 = p2.atoms
        self.assertEqual(atoms.cell, a2.cell)
        self.assertEqual(atoms.regions['left'].cell, a2.regions['left'].cell)
        self.assertEqual(len(atoms.regions['left']), len(a2.regions['left']))
        self.assertEqual(atoms.regions['right'].cell, a2.regions['right'].cell)
        self.assertEqual(len(atoms.regions['right']), len(a2.regions['right']))
        self.assertEqual(atoms.positions, a2.positions)
        self.assertEqual(str(atoms.symbols), str(a2.symbols))
