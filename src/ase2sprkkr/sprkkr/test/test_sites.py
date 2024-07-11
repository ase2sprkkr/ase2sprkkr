import numpy as np
if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..sites import SiteType
    from ..radial_meshes import ExponentialMesh
    from ..radial import RadialPotential, RadialCharge
    from ..sprkkr_atoms import SPRKKRAtoms


class TestSites(TestCase):

  def test_adding(self):
      atoms=SPRKKRAtoms(symbols='NK',pbc=True)
      atoms.cell = [[1,0,0], [0,1,1],[0,0,1]]
      assert atoms.sites[0].site_type != atoms.sites[1].site_type
      out = atoms + atoms
      assert len(out) == 4
      assert id(out.sites[0]) != id(atoms.sites[0])
      assert id(out.sites[0].site_type) != id(atoms.sites[0].site_type)
      half = out[:2]
      assert id(out.sites[0]) != id(half.sites[0])
      assert id(out.sites[0].site_type) != id(half.sites[0].site_type)

  def test_potential(self):
      a=SPRKKRAtoms('NaCl')
      m1 = ExponentialMesh(1.,0.005,None,None,200, None)
      m2 = ExponentialMesh(0.9,0.005,None,None,200, None)
      s0 = a.sites[0]
      s1 = a.sites[1]

      s0.potential = RadialPotential(m1.coors, m1, 2)
      s1.potential = s0.potential.copy()
      self.assertEqual(s0.potential(2.), 2.)
      self.assertEqual(s0.potential(1.), 1.)
      self.assertEqual(s0.potential(1.5), 1.5)
      self.assertEqual(s1.potential(2.), 2.)

      s1.potential = s0.potential.copy()
      self.assertEqual(s0.mesh, s1.mesh)

      s0.potential = p2 = RadialPotential(m2.coors * 2, m2, 2)
      self.assertEqual(s0.potential(2.), 4.)
      self.assertNotEqual(s0.mesh, s1.mesh)

      s0.mesh = m1
      self.assertEqual(s0.potential(2.), 4.)
      self.assertEqual(s0.potential(1.), 2.)
      self.assertEqual(s0.potential(1.5), 3.)

      self.assertEqual(s0.potential(2.), 4.)
      self.assertEqual(s0.mesh, s1.mesh)
      self.assertEqual(s0.potential.mesh, s1.mesh)
      s0.charge = m1.coors * 3
      self.assertEqual(s0.mesh, s1.mesh)
      self.assertEqual(s0.charge(2.), 6.)
      self.assertEqual(s0.charge(1.), 3.)
      s0.potential = p2
      self.assertEqual(s0.mesh, m2)
      self.assertEqual(s0.potential.mesh, m2)
      self.assertEqual(s0.charge.mesh, m2)
      self.assertEqual(s0.charge(2.), 6.)
      self.assertEqual(s0.charge(1.), 3.)
      s0.charge = RadialCharge(s1.potential._value, m1)
      self.assertEqual(s0.mesh, m1)
      self.assertEqual(s0.potential.mesh, m1)
      self.assertEqual(s0.charge.mesh, m1)

      s0.potential = [ m1.coors, m1.coors * 2 ]
      self.assertEqual(s0.potential.interpolate(1.5), np.asarray([1.5,3]))

      # test 1/r interpolation
      s0.potential = [ 1 / m1.coors, 1 / m1.coors * 2 ]
      self.assertEqual(s0.potential.interpolate(1.5), np.asarray([1 / 1.5,1 / 1.5 * 2]))

  def test_vacuum(self):
      a=SPRKKRAtoms('NaCl')
      site=SiteType(a, { 'Na' : 1.0 })
      self.assertFalse(site.is_vacuum())
      site=SiteType(a, { 'Na' : 0.5, 'Cl': 0.5 })
      self.assertFalse(site.is_vacuum())
      site=SiteType(a, { 'Na' : 0.5, 'Vc': 0.5 })
      self.assertFalse(site.is_vacuum())
      site=SiteType(a, { 'Vc': 1.0 })
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
      self.assertAlmostEqual(0.3 * 5. / 4., site.occupation[1])
      self.assertAlmostEqual(0.5 * 5. / 4., site.occupation[2])
      self.assertEqual(3, len(site.occupation))
      site.occupation.clean()
      self.assertEqual(2, len(site.occupation))
      self.assertAlmostEqual(0.3 * 5. / 4., site.occupation[0])
      self.assertAlmostEqual(0.5 * 5. / 4., site.occupation[1])
