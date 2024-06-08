import ase
import ase.build
import numpy as np
from ase2sprkkr import SPRKKR
import os

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from .. import symmetry as sy   # NOQA E402
from .. import spheres as sph   # NOQA E402


class TestSpheres(TestCase):

  def test(self):
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

      o = sph.empty_spheres(a2, point_symmetry =out)
      self.assertEqual(o.positions, np.asarray([[1.27279221,1.27279221,0.],[0,0,1.8]]))
      self.assertEqual(o.radii, np.asarray([1.44985909,1.44985909]))
      cu=ase.build.bulk('Cu')
      if os.environ.get('DO_NOT_RUN_SPRKKR', '') == '':
        out = SPRKKR().calculate(cu, empty_spheres=True, options={'niter': 1, 'ne' : 20, 'nktab' : 5 }, print_output=False)
        self.assertEqual(len(cu), 4)
