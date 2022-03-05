if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase import Atoms
from ase.build import bulk
from .. sprkkr_atoms import SPRKKRAtoms

class SPRKKRAtomsTest(TestCase):

 def xtest_atoms(self):
     a = Atoms('NaCl')
     a.set_positions([[0,0,0],[0,1,0]])
     a.info['occupancy'] = { 1: {'Cl' : 0.4, 'I' : 0.6 } }
     self.assertEqual(str(a.symbols), 'NaCl')
     SPRKKRAtoms.promote_ase_atoms(a)
     self.assertEqual(len(a.sites), 2)
     self.assertEqual(len(a.sites[0].occupation), 1)
     self.assertEqual(len(a.sites[1].occupation), 2)
     self.assertEqual(a.sites[1].occupation['Cl'], 0.4)
     self.assertEqual(str(a.symbols), 'NaI')
     self.assertEqual(a.sites[0].occupation.as_dict, {'Na' :  1})
     self.assertEqual(a.sites[1].occupation.as_dict, {'Cl' :  0.4, 'I': 0.6})
     self.assertEqual(a.info['occupancy'][1], {'Cl' :  0.4, 'I': 0.6})
     a.sites[1] = a.sites[0]
     a.compute_sites_symmetry()
     a.sites[1].occupation = 'Cl'
     self.assertEqual(str(a.symbols), 'NaCl')


 def test_symmetry(self):
     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     SPRKKRAtoms.promote_ase_atoms(atoms)
     self.assertTrue(atoms.sites[1] == atoms.sites[3])
     atoms.symmetry = False
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     atoms.symmetry = True
     self.assertTrue(atoms.sites[1] == atoms.sites[3])
     atoms.sites[3] = atoms.sites[3].copy()
     #No effect
     atoms.symmetry = True
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     atoms.compute_sites_symmetry()
     self.assertTrue(atoms.sites[1] == atoms.sites[3])
     atoms.symmetry = False
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     atoms.sites[3] = atoms.sites[1]
     #No effect
     atoms.symmetry = False
     self.assertTrue(atoms.sites[1] == atoms.sites[3])
     atoms.cancel_sites_symmetry()
     self.assertFalse(atoms.sites[1] == atoms.sites[3])

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     SPRKKRAtoms.promote_ase_atoms(atoms, symmetry=False)
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     #default None for symmetry do not change the already
     #initialized atoms object
     SPRKKRAtoms.promote_ase_atoms(atoms)
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     SPRKKRAtoms.promote_ase_atoms(atoms, symmetry=True)
     self.assertTrue(atoms.sites[1] == atoms.sites[3])


 def test_occupancy(self):
     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     SPRKKRAtoms.promote_ase_atoms(atoms)
     self.assertTrue(atoms.sites[1] == atoms.sites[3])
     atoms.cancel_sites_symmetry()
     self.assertFalse(atoms.sites[1] == atoms.sites[3])

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     atoms.info["occupancy"] = {
         0: { 'Li' : 1 },
         1: { 'Cl' : 1 },
         2: { 'Li' : 1 },
         3: { 'Cl' : 1 },
     }
     SPRKKRAtoms.promote_ase_atoms(atoms)
     self.assertTrue(atoms.sites[1] == atoms.sites[3])

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     atoms.info["occupancy"] = {
         0: { 'Li' : 1 },
         1: { 'Cl' : 0.5, 'I' :0.5 },
         2: { 'Li' : 1 },
         3: { 'Cl' : 1 },
     }
     SPRKKRAtoms.promote_ase_atoms(atoms)
     self.assertFalse(atoms.sites[1] == atoms.sites[3])
     self.assertEqual({ 'Cl' : 0.5, 'I' :0.5 }, atoms.sites[1].occupation.as_dict)
