if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..potentials import Potential
from ..potential_sections import PotentialSection
from ..potential_definitions import PotentialDefinition
from ase.spacegroup import crystal
import io
from ...common.grammar import generate_grammar

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

  def test_potential(self):
    a = 5.64
    nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
                   cellpar=[a, a, a, 90, 90, 90])

    sio = io.StringIO()
    p1 = Potential.from_atoms(nacl)
    p1.atoms.symbols[0] = 'K'
    p1.save_to_file(sio)
    sio.seek(0)
    p2 = Potential.from_file(sio)

    self.assertIsNotNone(p2.atoms)
    a1 = p1.atoms
    a2 = p2.atoms
    self.assertEqual(a2.cell * 1, a1.cell * 1)
    self.assertEqual(a2.positions, a1.positions)
    self.assertEqual(str(a2.symbols), str(a1.symbols))

  def test_examples(self):
    path = os.path.join(os.path.dirname(__file__), '../examples')
    i = 0

    #import cProfile, pstats, io
    #from pstats import Stats
    #with cProfile.Profile() as pr

    for x in os.listdir(path):
      if not x.endswith('.pot'): continue
      i+=1
      pot = Potential.from_file(os.path.join(path, x))
      self.assertTrue(isinstance(pot, Potential))
      self.assertTrue(pot.atoms.positions.shape[0] > 0)

    #Stats(pr).sort_stats('cumtime').print_stats(0.05)

    self.assertTrue(i >= 2)

  def test_reset(self):
    pot = Potential.from_file(os.path.join(os.path.dirname(__file__), '../examples/fp_new.pot'))
    self.assertNotEquals(pot.SCF_INFO.SCFSTATUS(), 'START')
    pot.reset()
    self.assertRaises(AttributeError, lambda: pot.CHARGE)
    self.assertRaises(AttributeError, lambda: pot.POTENTIAL)
    self.assertTrue(isinstance(pot.LATTICE, PotentialSection))
    self.assertEquals(pot.SCF_INFO.SCFSTATUS(), 'START')
