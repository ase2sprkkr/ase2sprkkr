if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase.build import bulk
from ..calculator import SprKkr
import os
import sys

class CalculatorTest(TestCase):

  def test_calculator(self):
     print_output = '-v' in sys.argv or '--verbose' in sys.argv

     #out = SprKkr.Potential.from_file(os.path.dirname(__file__) + '/output_test_calc.pot_new')

     atoms = bulk('Li')
     calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(options = {'NITER' : 2 }, print_output=print_output)
     self.assertEqual(2, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)

     #read again the output from a file - the results should be the same
     SprKkr.Task.create('scf').read_output_from_file(os.path.join(os.path.dirname(__file__), 'output_test_calc.out'))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(2, len(out.iterations))

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(options = {'NITER' : 1, 'NE' : 12 }, print_output=print_output)
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(3, len(out.iterations[-1]['atoms']))
     for i in out.iterations[-1]['atoms']:
        self.assertEqual(5, len(i['orbitals']))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))

     SprKkr.Task.create('scf').read_output_from_file(os.path.join(os.path.dirname(__file__), 'output_test_calc.out'))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))

     atoms = bulk('Li')
     calculator = SprKkr(atoms = atoms, mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     task = SprKkr.Task.create('scf')
     task.SCF.NITER = 1
     out = calculator.calculate(task = task, print_output=print_output)
     self.assertEqual(out.atoms, out.potential.atoms)
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)

     calculator = SprKkr(mpi=False, directory = os.path.dirname(__file__), input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(potential = out.potential, task = task)
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))
