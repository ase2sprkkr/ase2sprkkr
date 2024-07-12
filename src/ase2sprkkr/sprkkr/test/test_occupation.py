import os

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


class test(TestCase):

    def test(self):
        from ase import io
        from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms
        ciffile='test.cif'
        atoms = io.read(os.path.join(
            os.path.dirname(__file__),
            ciffile
        ))
        SPRKKRAtoms.promote_ase_atoms(atoms)
        self.assertEqual(atoms.sites[0].occupation['Hf'], 0.7025)
        self.assertEqual(atoms.sites[0].occupation['Ti'], 0.2975)
        self.assertEqual(atoms.sites[7].occupation['H'], 0.233)
        self.assertEqual(atoms.sites[7].occupation['H'], 0.233)
