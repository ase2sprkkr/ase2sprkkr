import os
import pytest
from ase.build import bulk

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

class Test(TestCase):
    def test_occ(self):
        from ase import io
        from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms

        ciffile = "test.cif"
        atoms = io.read(os.path.join(os.path.dirname(__file__), ciffile))
        SPRKKRAtoms.promote_ase_atoms(atoms)
        self.assertEqual(atoms.sites[0].occupation["Hf"], 0.7025)
        self.assertEqual(atoms.sites[0].occupation["Ti"], 0.2975)
        self.assertEqual(atoms.sites[-1].occupation["H"], 0.233)

        atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
        SPRKKRAtoms.promote_ase_atoms(atoms)
        with pytest.raises(RuntimeError):
           atoms.sites[0].occupation['Fe'] = 0.5
        atoms.sites[0].site_type.occupation['Fe'] = 0.5
        atoms.sites[0].break_symmetry()
        atoms.sites[0].occupation['Fe'] = 0.5
