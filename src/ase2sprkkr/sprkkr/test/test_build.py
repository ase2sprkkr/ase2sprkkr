if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..build import semiinfinite_system  # NOQA E402
from ase.build import bulk               # NOQA E402


class TestBuild(TestCase):

    def test(self):
        atoms=bulk('LiCl', 'rocksalt', a=5.64)
        si=semiinfinite_system(atoms,(4.2))
        self.assertEqual(len(si), 22)
        self.assertEqual(sum(si.symbols == 'Li'), 6)
        self.assertEqual(sum(si.symbols == 'Cl'), 5)
        self.assertEqual(sum(si.symbols == 'X'), 11)
        self.assertEqual(si.positions, (atoms * (1,1,11)).positions)
