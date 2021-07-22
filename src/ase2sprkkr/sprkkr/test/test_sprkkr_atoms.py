if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase import Atoms
from .. sprkkr_atoms import SPRKKRAtoms

class SPRKKRAtomsTest(TestCase):

 def test_calculator(self):
     a = Atoms('NaCl')
     a.set_positions([[0,0,0],[0,1,0]])
     a.info['occupancy'] = { 1: {'Cl' : 0.4, 'I' : 0.6 } }
     self.assertEquals(str(a.symbols), 'NaCl')
     SPRKKRAtoms.promote_ase_atoms(a)
     self.assertEquals(len(a.sites), 2)
     self.assertEquals(len(a.sites[0].occupation), 1)
     self.assertEquals(len(a.sites[1].occupation), 2)
     self.assertEquals(a.sites[1].occupation['Cl'], 0.4)
     self.assertEquals(str(a.symbols), 'NaI')
     self.assertEquals(a.sites[0].occupation.as_dict, {'Na' :  1})
     self.assertEquals(a.sites[1].occupation.as_dict, {'Cl' :  0.4, 'I': 0.6})
     self.assertEquals(a.info['occupancy'][1], {'Cl' :  0.4, 'I': 0.6})
     a.sites[1] = a.sites[0]
     a.compute_sites_symmetry()
     a.sites[1].occupation = 'Cl'
     self.assertEquals(str(a.symbols), 'NaCl')
