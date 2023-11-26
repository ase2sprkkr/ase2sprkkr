if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from os import listdir
from os.path import dirname, join as path_join, isfile
import os.path
from ..output_files import OutputFile

import scipy.constants
from ase.units import Rydberg as Ry
from ...common.configuration_containers import DisabledAttributeError
import tempfile

class TestOutput(TestCase):

  def test_output(self):
    dire = path_join(dirname(dirname(__file__)), 'output_files', 'examples')
    for i in listdir(dire):
        fname = path_join(dire, i)
        if isfile(fname):
            ext=fname.rsplit('.',1)[1]
            out=OutputFile.from_file(fname, unknown=False)
            if ext=='spc':
               self.assertEqual('ARPES', out.KEYWORD())
               self.assertEqual((200,160), out.ENERGY().shape)
               o2 = out + out
               for i in 'ENERGY', 'THETA', 'K', 'DETERMINANT':
                   self.assertEqual(out[i](), o2[i]())
               for i in 'TOTAL', 'POLARIZATION', 'UP', 'DOWN':
                   self.assertEqual(2*out[i](), o2[i]())

            elif ext=='dos':
               self.assertEqual(out.n_orbitals(1), 3)
               self.assertEqual(out.n_spins(), 2)
               self.assertEqual((2,3,1200), out.dos_for_site_type('Ta').shape)
               self.assertEqual(out.DOS['Ta'][5] / Ry, out.dos_for_site_type('Ta',1,2)[:])

            elif ext=='bsf':
               self.assertEqual(out.I().shape, (out.NQ_EFF(), out.NE(), out.NK()))
               if out.KEYWORD() == 'BSF':
                  self.assertEqual(out.I_UP().shape, (out.NQ_EFF(), out.NE(), out.NK()))
                  self.assertRaises(DisabledAttributeError, lambda: out.I_X)
               else:
                  self.assertEqual(out.I_X().shape, (out.NQ_EFF(), out.NE(), out.NK()))
                  self.assertRaises(DisabledAttributeError, lambda: out.I_UP)

               if out.MODE() == 'EK-REL':
                  self.assertEqual(len(out.K()), out.NK())
                  self.assertEqual(len(out.E()), out.NE())
               else:
                  self.assertWarns(RuntimeWarning, out.NK1())
                  self.assertEqual(out.NK1()[0], 0)
                  out.VECK_START=[1,1,1]
                  self.assertEqual(out.NK1()[0], np.sqrt(3))
                  self.assertEqual(len(out.K1()), out.NK1())
                  self.assertEqual(len(out.K2()), out.NK2())

            with tempfile.NamedTemporaryFile() as name:
                  out.plot(filename=name)
