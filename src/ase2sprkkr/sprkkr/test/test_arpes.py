import os
import sys
import numpy as np
from ase.atoms import Atoms


if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


from ..calculator import SPRKKR   # NOQA
from ...potentials.potentials import Potential  # NOQA


class TestCalculator(TestCase):

 print_output = '-v' in sys.argv or '--verbose' in sys.argv
 dirname = os.path.dirname(__file__)
 _calc_args = dict(
     directory = dirname, input_file = 'output_test_calc.inp',
     output_file = 'output_test_calc.out',
     potential_file ='output_test_calc.pot',
     print_output=print_output,
     mpi = 'auto',
     empty_spheres=False,
 )

 @classmethod
 def calc_args(cls, **kwargs):
      kwargs.update(cls._calc_args)
      return kwargs

 def run_sprkkr(self):
      return os.environ.get('DO_NOT_RUN_SPRKKR', '') == ''

 def test_run(self):
      if not self.run_sprkkr():
          return

      a=Atoms(symbols='Cu',cell=np.array([(1.,0,0),(0,1,0),(0,0,1)]))
      a.pbc=True
      xx=SPRKKR(atoms=a, **self.calc_args())
      xx.input_parameters.SCF.NITER=1
      xx.set('NE', 10)
      xx.set('NKTAB', 10)
      out = xx.calculate()
      calc = out.calculator
      calc.change_task('arpes')
      out2 = out.calculator.calculate()
      self.assertTrue(isinstance(out2.spc.ENERGY(), np.ndarray))
