import os
import io
from ase.spacegroup import crystal
from datetime import datetime
import numpy as np

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..potentials import Potential
    from ..potential_sections import PotentialSection
    from ...sprkkr.calculator import SPRKKR


class TestPotential(TestCase):

  def test_grammar(self):
    g=Potential.potential_definition.grammar()

    def check(x, prefix=''):
       if '\n' in x.whiteChars:
          raise Exception(f'The potential grammar has newline in whitespace, which is illegal:\n {x}: {prefix}:{x.whiteChars}')
       if hasattr(x, 'expr'):
          check(x.expr, prefix + '.expr')
       if hasattr(x, 'exprs'):
          for i,e in enumerate(x.exprs):
              check(e, prefix + f'.exprs[{i}]')
    check(g)

  def test_potential(self, temporary_dir):
    a = 5.64
    nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
                   cellpar=[a, a, a, 90, 90, 90])

    sio = io.StringIO()
    p1 = Potential.from_atoms(nacl)
    p1.atoms.symbols[0] = 'K'
    p1.save_to_file(sio)
    sio.seek(0)
    p2 = Potential.from_file(sio)

    assert p2.atoms is not None
    a1 = p1.atoms
    a2 = p2.atoms
    self.assertEqual(a2.cell * 1, a1.cell * 1)
    self.assertEqual(a2.positions, a1.positions)
    self.assertEqual(str(a2.symbols), str(a1.symbols))

  def test_2D(self, temporary_dir):
    path = os.path.join(os.path.dirname(__file__), '..','examples','GeTe.pot')
    p=Potential.from_file(path)
    if os.environ.get('DO_NOT_RUN_SPRKKR', '') == '':
        SPRKKR().calculate(potential=p, **self.calc_args(options={'NITER':1, 'NKTAB':5, 'NE':5}))
    else:
        SPRKKR().save_input(potential=p, options={'NITER':1}, **self.calc_args())

  def test_potential_data(self, temporary_dir):
    path = os.path.join(os.path.dirname(__file__), '..','examples','FePt.new.pot')
    p=Potential.from_file(path)
    self.assertTrue(p.atoms.sites[0].potential is not None)
    self.assertTrue(p.atoms.sites[0].charge is not None)
    assert p.atoms.sites[-1].potential is not None
    assert p.atoms.sites[-1].charge is not None
    p.atoms.sites[1].potential.bt = -p.atoms.sites[1].potential.bt
    s=io.StringIO()
    p.save_to_file(s)
    s.seek(0)
    pp = Potential.from_file(s)
    p=Potential.from_file(path)
    self.assertEqual(pp.atoms.sites[1].potential.bt, -p.atoms.sites[1].potential.bt)
    self.assertFalse(np.allclose(pp.atoms.sites[1].potential.bt, p.atoms.sites[1].potential.bt))

    if os.environ.get('DO_NOT_RUN_SPRKKR', '') == '':
        SPRKKR().calculate(potential=p, **self.calc_args(options={'NITER':1,'NKTAB':5, 'NE':5}))

  def test_examples(self):
    path = os.path.join(os.path.dirname(__file__), '..','examples')
    i = 0

    # import cProfile, pstats, io
    # from pstats import Stats
    # with cProfile.Profile() as pr

    for x in os.listdir(path):
      if not x.endswith('.pot'):
          continue
      i+=1
      pot = Potential.from_file(os.path.join(path, x))
      self.assertTrue(isinstance(pot, Potential))
      self.assertTrue(pot.atoms.positions.shape[0] > 0)

    # Stats(pr).sort_stats('cumtime').print_stats(0.05)

    self.assertTrue(i >= 2)

  def test_sequence(self):
    p = Potential()
    self.assertEqual(7, p.HEADER.FORMAT.VERSION)
    p.HEADER.FORMAT.VERSION=5
    self.assertEqual(5, p.HEADER.FORMAT.VERSION)
    p.HEADER.FORMAT[0] = 8
    self.assertEqual(8, p.HEADER.FORMAT.VERSION)
    self.assertEqual(datetime, p.HEADER.FORMAT.DATE.__class__)
    self.assertEqual(datetime, p.HEADER.FORMAT[1].__class__)

  def test_reset(self):
    pot = Potential.from_file(os.path.join(os.path.dirname(__file__), '../examples/fp_new.pot'))
    pot.atoms.sites[0].mesh.r1 = 1.
    pot.MESH_INFORMATION.DATA()[0][0]=2.
    self.assertNotEqual(pot.SCF_INFO.SCFSTATUS(), 'START')
    # the reset of the potential is just reseting all non-mandatory
    # properties. The mandatory are currently not set to defaults!
    pot.reset()
    self.assertEqual(2., pot.atoms.sites[0].mesh.r1)
    self.assertEqual(2., pot.MESH_INFORMATION.DATA()[0][0])
    self.assertEqual(len(pot.CHARGE), 6)
    self.assertEqual(len(pot.POTENTIAL), 6)
    self.assertTrue(isinstance(pot.LATTICE, PotentialSection))
    self.assertEqual(pot.SCF_INFO.SCFSTATUS(), 'START')
    # a more hard version - it resets all the informations
    pot.atoms.reset_sprkkr_potential()
    self.assertEqual(pot.MESH_INFORMATION.DATA()[0][0], pot.atoms.sites[0].mesh.r1)
    self.assertNotEqual(1., pot.atoms.sites[0].mesh.r1)
    self.assertNotEqual(2., pot.MESH_INFORMATION.DATA()[0][0])
