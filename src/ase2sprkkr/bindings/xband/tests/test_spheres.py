import ase
import ase.build
import numpy as np
from ase2sprkkr import SPRKKR
from ase.units import Bohr
import os

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from .. import symmetry as sy   # NOQA E402
from .. import spheres as sph   # NOQA E402
from ...empty_spheres import add_empty_spheres, empty_spheres # NOQA E402
from .... import Potential       # NOQA E402


class TestSpheres(TestCase):

  def test_xband(self):
      dirr = os.path.dirname(__file__)
      pot = Potential.from_file(os.path.join(dirr, 'MnTi3.pot'))
      mnti= pot.atoms
      sym = sy.find_symmetry(mnti)
      sym2 = sy.find_symmetry(mnti, use_spacegroup=False)
      self.assertEqual(sym, sym2)

      pot = Potential.from_file(os.path.join(dirr, 'Cu.pot'))
      full = pot.atoms
      assert full.sites[0].site_type is full.sites[1].site_type
      cu = pot.atoms[:2]
      assert cu.sites[0].site_type is cu.sites[1].site_type
      assert cu.spacegroup_info.number() == 227
      cu.sites[0].break_symmetry()
      assert cu.sites[0].site_type is not cu.sites[1].site_type
      cu.compute_sites_symmetry(consider_old=False)
      assert cu.sites[0].site_type is cu.sites[1].site_type

      sym = sy.find_symmetry(cu)
      should = np.array([[
        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
        33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56],[
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],
        dtype=np.int32)
      self.assertEqual(sym, should)
      sym = sy.find_symmetry(cu, use_spacegroup=False)
      self.assertEqual(sym, should)
      add_empty_spheres(cu, method='xband')
      sort = lambda x: np.asarray(sorted(tuple(i) for i in x))
      self.assertEqual(sort(cu.get_scaled_positions()), sort(full.get_scaled_positions()))

      pot = Potential.from_file(os.path.join(dirr, 'V.pot'))
      v = pot.atoms
      sym = sy.find_symmetry(v)
      sym2 = sy.find_symmetry(v, subprocess=False, use_spacegroup=False)
      self.assertEqual(sym, sym2)
      self.assertEqual(len(empty_spheres(v, method='xband')), 0)

  def test(self, temporary_dir):
      a2 = ase.build.bulk('Cu', 'fcc', a=3.6, orthorhombic=True)
      sg = ase.spacegroup.get_spacegroup(a2)
      cp = a2.cell.cellpar()
      sy.find_symmetry_ex(sg.no, cp[:3],cp[3:], a2.cell[:], len(a2), np.ascontiguousarray(a2.get_scaled_positions()), len(a2), np.arange(len(a2), dtype=np.int32), a2.positions, False)
      out = sy.find_symmetry(a2, subprocess=False)
      self.assertEqual(out.shape, (2,48))
      out = sy.find_symmetry(a2, subprocess=True)
      self.assertEqual(out.shape, (2,48))
      out2 = sy.find_symmetry(a2, use_spacegroup=False, subprocess=False)
      self.assertEqual(out, out2)

      # these numbers should be returned, at least I hope so
      self.assertEqual(out[0,47], 56)
      self.assertEqual(out[1,47], 1)

      o = sph.empty_spheres(a2, point_symmetry=out, min_radius = 0.7)
      self.assertEqual(len(o), 0)
      o = sph.empty_spheres(a2, point_symmetry=out, min_radius = 0.5)
      self.assertEqual(len(o), 10)
      self.assertEqual(o.radii, np.asarray([0.995084964973998 * Bohr] * 10))
      cu=ase.build.bulk('Cu')
      if os.environ.get('DO_NOT_RUN_SPRKKR', '') == '':
        out = SPRKKR().calculate(cu, **self.calc_args(empty_spheres={'min_radius': 0.25}, options={'niter': 1, 'ne' : 20, 'nktab' : 5 }))
        self.assertEqual(len(cu), 4)
