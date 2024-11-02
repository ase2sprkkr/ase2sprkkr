from ase.build import bulk
import numpy as np

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ...sprkkr.sprkkr_atoms import SPRKKRAtoms, SpacegroupInfo  # NOQA: E402


class Test(TestCase):

  def test_spacegroup_info(self):

      def ssert(mapping, out):
           self.assertEqual(out, np.asarray(mapping))

      atoms = bulk('NaCl', "rocksalt", a=5.64)
      sgi = SpacegroupInfo(atoms)
      ssert([0,1], sgi.equivalent_sites)
      self.assertEqual(225, sgi.spacegroup_number())
      self.assertEqual(225, sgi.dataset.number)

      atoms = SPRKKRAtoms('ION')
      atoms.set_positions(np.arange(9).reshape((3,3)))
      sgi = SpacegroupInfo(atoms)
      assert sgi.dataset is False
      ssert([0,1,2], sgi.equivalent_sites)
