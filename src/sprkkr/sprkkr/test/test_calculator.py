if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase.build import bulk
from ..calculator import SprKkr
import os

class CalculatorTest(TestCase):

  def test_calculator(self):
      atoms = bulk('Na')
      #atoms.cell = atoms.cell / 2 #lwoer computational requirements
      calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
      calculator.calculate(options = {'NITER' : 1 }, print_output=True)
