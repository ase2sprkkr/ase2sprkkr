if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase.build import bulk
from ..calculator import SPRKKR
from ...potentials.potentials import Potential
import os
import sys
from pathlib import Path

class CalculatorTest(TestCase):

 def test_calculator(self):
     print_output = '-v' in sys.argv or '--verbose' in sys.argv
     dirname = os.path.dirname(__file__)
     here = lambda x: os.path.join(dirname, x)

     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot', print_output=print_output)
     calculator.input_parameters.set(NE=21111)
     self.assertEqual(calculator.input_parameters.get('NE'),21111)
     calculator.set(NE=31111)
     self.assertEqual(calculator.input_parameters.get('NE'),31111)
     self.assertEqual(calculator.input_parameters.get('ENERGY.NE'),31111)
     calculator.input_parameters.set({'ENERGY.NE':41111})
     self.assertEqual(calculator.input_parameters.get('NE'),41111)
     calculator.set({'ENERGY.NE':51111})
     self.assertEqual(calculator.input_parameters.get('NE'),51111)
     calculator.set('ENERGY.NE', 451111)
     self.assertEqual(calculator.input_parameters.get('NE'),451111)

     calculator.input_parameters.set(NE=11111)
     self.assertEqual(calculator.input_parameters.get('NE'),11111)
     self.assertEqual(calculator.get('NE'),11111)
     self.assertEqual(calculator.input_parameters.ENERGY.NE(),11111)
     self.assertEqual(calculator.input_parameters.TASK.TASK(), 'SCF')
     calculator.input_parameters = 'PHAGEN'
     self.assertEqual(calculator.input_parameters.TASK.TASK(), 'PHAGEN')
     self.assertNotEqual(calculator.input_parameters.get('NE'), 111111)

     calculator = SPRKKR(atoms = atoms, mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot', print_output=print_output)

     inp_file=here('output_test_calc.inp')
     pot_file=here('output_test_calc.pot')
     Path(inp_file).touch()
     Path(pot_file).touch()

     def assert_change(*args):
        for ext, res in zip(('inp','pot'), args):
          fname = here('output_test_calc.' + ext)
          ntime = os.path.getmtime(fname)
          changed = abs(ntime -  time[ext]) > 2.0
          self.assertEqual(res, changed)
          if changed:
            ntime -= 10
            time[ext] = ntime
            os.utime(fname, ( ntime, ntime ))

     time = { 'inp' : 0, 'pot' : 0 }
     assert_change(True, True)
     assert_change(False, False)
     calculator.save_input()
     assert_change(True, True)
     assert_change(False, False)
     calculator.save_input(input_parameters = inp_file)
     assert_change(True, True)
     calculator.save_input(input_parameters = inp_file, potential = False)
     assert_change(False, False)
     calculator.save_input(input_parameters = inp_file, options = {'NITER' : 2}, potential = False)
     assert_change(True, False)
     calculator.save_input(input_parameters = inp_file, potential = atoms.potential)
     assert_change(True, True)
     calculator.save_input(input_parameters = inp_file, potential = pot_file)
     assert_change(True, False)

 def test_run(self):
     ignore = os.environ.get('DO_NOT_RUN_SPRKKR', '') != ''
     if ignore:
        return

     print_output = '-v' in sys.argv or '--verbose' in sys.argv
     dirname = os.path.dirname(__file__)
     here = lambda x: os.path.join(dirname, x)

     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot', print_output=print_output)
     #use methods of atoms
     atoms.calc = calculator
     self.assertTrue(isinstance(atoms.get_potential_energy(), float))

     #calculator options
     out = calculator.calculate(options = {'NITER' : 2 }, print_output=print_output, mpi=4)
     self.assertEqual(2, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)
     self.assertTrue(isinstance(out.potential, Potential))
     self.assertTrue(isinstance(out.calculator.potential_object, Potential))

     #read again the output from a file - the results should be the same
     out = SPRKKR.InputParameters.create('scf').read_output_from_file(here('output_test_calc.out'))
     self.assertEqual(2, len(out.iterations))
     out.plot(filename = here('output_test_calc.png'))

     #use methods of atoms
     atoms.calc = calculator
     self.assertTrue(isinstance(atoms.get_potential_energy(), float))

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     calculator = SPRKKR(atoms = atoms, mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(options = {'NITER' : 1, 'NE' : 12 }, print_output=print_output)
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(3, len(out.iterations[-1]['atoms']))
     for i in out.iterations[-1]['atoms']:
        self.assertEqual(4, len(i['orbitals']))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))

     out = SPRKKR.InputParameters.from_file(here('output_test_calc.inp')).read_output_from_file(here('output_test_calc.out'))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))

     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     ips = SPRKKR.InputParameters.create('scf')
     ips.SCF.NITER = 1
     out = calculator.calculate(input_parameters=ips, print_output=print_output)
     self.assertEqual(out.atoms, out.potential.atoms)
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)

     calculator = SPRKKR(mpi=False, directory = dirname, input_file = 'output_test_calc.inp', output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot')
     out = calculator.calculate(potential=out.potential,input_parameters=ips,print_output=print_output)
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))

     calculator.calculate(input_parameters='PHAGEN', potential=out.potential_filename)
     out.calculator.calculate(input_parameters='PHAGEN')
