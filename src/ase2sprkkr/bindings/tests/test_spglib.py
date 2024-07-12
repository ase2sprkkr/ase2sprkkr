from ase.build import bulk
import numpy as np

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..spglib import possibly_equivalent_sites, SpacegroupInfo  # NOQA: E402
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms                  # NOQA: E402


class TestSpgilib(TestCase):

  def test_possibly_equivalent_sites(self):

       def ssert(mapping, out):
           self.assertTrue(out.is_equivalent_to(mapping))

       atoms = bulk('NaCl', "rocksalt", a=5.64)
       ssert([0,1], possibly_equivalent_sites(atoms))
       ssert([0,0], possibly_equivalent_sites(atoms, atomic_numbers=[4,4]))

       atoms = bulk('NaNa', "rocksalt", a=5.64)
       ssert([0,0], possibly_equivalent_sites(atoms))
       ssert([0,1], possibly_equivalent_sites(atoms, atomic_numbers=[4,7]))

       atoms = bulk('NaCl', "rocksalt", a=5.64)
       SPRKKRAtoms.promote_ase_atoms(atoms)
       atoms.sites
       ssert([0,1], possibly_equivalent_sites(atoms))
       atoms.sites[1].occupation = {'Na':1.0}
       ssert([0,1], possibly_equivalent_sites(atoms, consider_old=True))
       # repeat the test to be sure
       ssert([0,1], possibly_equivalent_sites(atoms, consider_old=True))
       ssert([0,0], possibly_equivalent_sites(atoms, consider_old=False))

  def test_spacegroup_info(self):

      def ssert(mapping, out):
           self.assertTrue(out.is_equivalent_to(mapping))

      atoms = bulk('NaCl', "rocksalt", a=5.64)
      sgi = SpacegroupInfo.from_atoms(atoms)
      ssert([0,1], sgi.equivalent_sites)
      self.assertEqual(225, sgi.spacegroup.no)
      self.assertEqual(225, sgi.number())
      self.assertEqual(225, sgi.dataset['number'])

      sgi = SpacegroupInfo(atoms, sgi.spacegroup)
      self.assertEqual(225, sgi.spacegroup.no)
      self.assertEqual(225, sgi.dataset['number'])
      ssert([0,1], sgi.equivalent_sites)

      atoms = SPRKKRAtoms('ION')
      atoms.set_positions(np.arange(9).reshape((3,3)))
      sgi = SpacegroupInfo.from_atoms(atoms)
      self.assertIsNone(sgi.spacegroup)
      self.assertIsNone(sgi.dataset)
      ssert([0,1,2], sgi.equivalent_sites)
