if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase.build import bulk
from ..calculator import SprKkr
import os
import sys


from ...calcio import PotFile

class CalculatorTest(TestCase):

  def test_calculator(self):
      #atoms = bulk('Li')
      #atoms.cell = atoms.cell / 2 #lwoer computational requirements
      #pf = PotFile(atoms=atoms,filename = 'output2_test_calc.pot', directory = os.path.dirname(__file__))
      #pf.write()
      #import subprocess
      #subprocess.run("kkrscf < output2_test_calc.inp", shell=True)
     print_output = '-v' in sys.argv or '--verbose' in sys.argv

     atoms = bulk('Li')
     calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(options = {'NITER' : 2 }, print_output=print_output)
     self.assertEquals(2, len(out))

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(options = {'NITER' : 1, 'NE' : 12 }, print_output=print_output)
     self.assertEquals(1, len(out))
     self.assertEquals(3, len(out[-1]['atoms']))
     for i in out[-1]['atoms']:
        self.assertEquals(5, len(i['orbitals']))


