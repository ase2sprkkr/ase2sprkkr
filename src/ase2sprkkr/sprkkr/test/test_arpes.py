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

 def test_run(self, temporary_dir):
      if not self.run_sprkkr():
          return

      a=Atoms(symbols='Cu',cell=np.array([(1.,0,0),(0,1,0),(0,0,1)]))
      a.pbc=True
      xx=SPRKKR(atoms=a, **self.calc_args())
      xx.input_parameters.SCF.NITER=1
      xx.set('NE', 10)
      xx.set('NKTAB', 10)
      out = xx.calculate(**self.calc_args())
      calc = out.calculator
      calc.change_task('arpes')
      out2 = out.calculator.calculate(**self.calc_args())
      self.assertTrue(isinstance(out2.spc.ENERGY(), np.ndarray))
