from ase.build import bulk
import os
import re
from pathlib import Path
from ase import Atoms

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


if True:
    from ..calculator import SPRKKR
    from ...potentials.potentials import Potential
    from ..build import semiinfinite_system


def _fast_atoms(b, jrws=20,r1=2e-6):
     for i in b.sites:
         i.mesh.jrws = jrws
         i.mesh.r1 =r1


class TestCalculator(TestCase):

 def test_2D(self, temporary_dir):
     a=Atoms(symbols="C", positions=[[0,0,0]], cell=[[1,0,0],[0,1,0], [0,0,1]], pbc=[1,1,1])
     b=semiinfinite_system(a, repeat=2)
     cal=SPRKKR(atoms=b, **self._calc_args)
     _fast_atoms(b)
     cal.input_parameters.set_from_atoms(b)
     self.assertTrue(bool(re.search('NKTAB3D=5', cal.input_parameters.to_string())))
     self.assertFalse(bool(re.match('NKTAB=', cal.input_parameters.to_string())))
     if not self.run_sprkkr():
         return
     out=SPRKKR()
     out=out.calculate(b, **self._calc_args(options={'EMIN': 8., 'NITER':2, 'NKTAB3D':2}))
     self.assertTrue(bool(re.search('NKTAB3D=', out.input_parameters.to_string())))
     self.assertFalse(bool(re.match('NKTAB=', out.input_parameters.to_string())))

 def test_calculator(self, temporary_dir):
     here = lambda x: os.path.join(self.dirname, x)

     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, **self.calc_args())
     _fast_atoms(atoms)
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

     calculator = SPRKKR(atoms = atoms, **self.calc_args())

     inp_file=here('output_test_calc.inp')
     pot_file=here('output_test_calc.pot')
     Path(inp_file).touch()
     Path(pot_file).touch()

     def assert_change(*args):
        for ext, res in zip(('inp','pot'), args):
          fname = here('output_test_calc.' + ext)
          ntime = os.path.getmtime(fname)
          changed = abs(ntime - time[ext]) > 2.0
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

 def test_run(self, temporary_dir):
     if not self.run_sprkkr():
         return

     here = lambda x: os.path.join(self.dirname, x)

     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, **self.calc_args(mpi=False))
     _fast_atoms(atoms)
     # use methods of atoms
     atoms.calc = calculator
     calculator.input_parameters.find('NITER').set(2)
     calculator.input_parameters.find('NKTAB').set(50)
     self.assertTrue(isinstance(atoms.get_potential_energy(), float))
     calculator.input_parameters.find('NITER').set(100)

     # calculator options
     out = calculator.calculate(options={'NITER' : 2 }, mpi=4)
     self.assertEqual(2, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)
     self.assertTrue(isinstance(out.potential, Potential))
     self.assertTrue(isinstance(out.calculator.potential_object, Potential))

     # read again the output from a file - the results should be the same
     out2 = SPRKKR.InputParameters.create('scf').read_output_from_file(out.output_file)
     self.assertEqual(out.output_file, out2.output_file)
     self.assertEqual(2, len(out2.iterations ))
     out.plot(filename = here('output_test_calc.png'))

     # use methods of atoms
     atoms.calc = calculator
     self.assertTrue(isinstance(atoms.get_potential_energy(), float))

     atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
     calculator = SPRKKR(atoms = atoms, **self.calc_args())
     # _fast_atoms(atoms, jrws=50, r1=1e-6)
     out = calculator.calculate(options = {'EMIN': -1., 'NITER' : 1, 'NKTAB': 40, 'NE' : 20})
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(3, len(out.iterations[-1]['atomic_types']))
     for i in out.iterations[-1]['atomic_types'].values():
        self.assertEqual(4, len(i['orbitals']))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))

     out = SPRKKR.InputParameters.from_file(out.input_file).read_output_from_file(out.output_file)
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))

 def test_phagen(self, temporary_dir):
     if not self.run_sprkkr():
         return
     atoms = bulk('Li')
     calculator = SPRKKR(atoms = atoms, **self.calc_args())
     _fast_atoms(atoms)
     ips = SPRKKR.InputParameters.create('scf')
     ips.SCF.NITER = 1
     ips.TAU.NKTAB=20
     out = calculator.calculate(input_parameters=ips)
     self.assertEqual(out.atoms, out.potential.atoms)
     self.assertEqual(1, len(out.iterations))
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertFalse(out.converged)

     calculator = SPRKKR(**self.calc_args())
     out = calculator.calculate(potential=out.potential,input_parameters=ips)
     self.assertEqual(str(atoms.symbols), str(out.atoms.symbols))
     self.assertEqual(1, len(out.iterations))

     calculator.calculate(input_parameters='PHAGEN', potential=out.potential_filename, options={'NKTAB' : 30})
     out.calculator.calculate(input_parameters='PHAGEN', options={'NKTAB' : 10})
