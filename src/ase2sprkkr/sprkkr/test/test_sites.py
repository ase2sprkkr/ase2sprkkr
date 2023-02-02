if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..sites import Site
from ..sprkkr_atoms import SPRKKRAtoms

class SitesTest(TestCase):

  def test_vacuum(self):
      a=SPRKKRAtoms('NaCl')
      site=Site(a, { 'Na' : 1.0 })
      self.assertFalse(site.is_vacuum())
      site=Site(a, { 'Na' : 0.5, 'Cl': 0.5 })
      self.assertFalse(site.is_vacuum())
      site=Site(a, { 'Na' : 0.5, 'Vc': 0.5 })
      self.assertFalse(site.is_vacuum())
      site=Site(a, { 'Vc': 1.0 })
      self.assertTrue(site.is_vacuum())

      a=SPRKKRAtoms('NaX')
      self.assertFalse(a.sites[0].is_vacuum())
      self.assertTrue(a.sites[1].is_vacuum())
      a.sites[1].atomic_types[0]='Li'
      self.assertFalse(a.sites[1].is_vacuum())

  def test_occupancy(self):
      a=SPRKKRAtoms('Na')
      site = a.sites[0]
      site.occupation = { 'Na' : 0.7, 'Cl' : None }
      self.assertAlmostEqual(0.3, site.occupation[1])
      self.assertAlmostEqual(0.3, site.occupation['Cl'])
      site.occupation['Cl'] = 0.6
      self.assertEqual('Cl', a.symbols)
      self.assertAlmostEqual(0.4, site.occupation[0])
      site.occupation['N'] = 0.5

      self.assertAlmostEqual(0.2, site.occupation[0])
      self.assertAlmostEqual(0.3, site.occupation[1])
      self.assertAlmostEqual(0.5, site.occupation[2])
      self.assertAlmostEqual(0.5, site.occupation['N'])
      self.assertEqual('N', a.symbols)
      site.occupation[0] = 0.
      self.assertEqual(0., site.occupation[0])
      self.assertAlmostEqual(0.3*5./4., site.occupation[1])
      self.assertAlmostEqual(0.5*5./4., site.occupation[2])
      self.assertEqual(3, len(site.occupation))
      site.occupation.clean()
      self.assertEqual(2, len(site.occupation))
      self.assertAlmostEqual(0.3*5./4., site.occupation[0])
      self.assertAlmostEqual(0.5*5./4., site.occupation[1])
