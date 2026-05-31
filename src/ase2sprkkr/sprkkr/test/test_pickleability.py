import pickle

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..calculator import SPRKKR
    from ..sprkkr_atoms import SPRKKRAtoms


class TestSPRKKRPickle(TestCase):

    def test_calculator_pickles(self):
        pickle.dumps(SPRKKR())

    def test_atoms_roundtrip(self):
        atoms = SPRKKRAtoms('H')
        atoms.calc = SPRKKR()

        restored = pickle.loads(pickle.dumps(atoms))

        self.assertEqual(restored.symbols[0], 'H')
        assert isinstance(restored.calc, SPRKKR)
