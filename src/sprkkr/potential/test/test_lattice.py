from .init_tests import TestCase, patch_package
__package__ = patch_package(__package__)

import os
from ..lattice_section import LatticeSection, LatticeSectionDefinition
from ..potential import AtomsWrapper
from ..potential_definitions import PotentialDefinition
from ase.build import bulk
import io


class TestLattice(TestCase):

  def test_lattice(self):
    pot = AtomsWrapper()
    df=LatticeSectionDefinition('LATTICE')
    sec=LatticeSection(df, pot)
    pd=PotentialDefinition('pot', [df])

    for e in ['Co', 'Ni', 'U', 'Ge', 'As', 'Na', 'In', 'Pr']:
     #for e in ['As']:
       atoms = bulk(e)
       pot.atoms = atoms
       cell = atoms.cell
       sec.cell = cell
       s = io.StringIO()
       sec.save_to_file(s)
       self.assertTrue(sec.cell is cell)
       s.seek(0)
       out = pd.read_from_file(s)
       import numpy as np
       from ...common.misc import OrderedDict

       a=out.LATTICE.to_dict()
       b=sec.to_dict()

       self.assertEqual(out.LATTICE.to_dict(), sec.to_dict())
       self.assertTrue(sec.cell is cell)
       self.assertIsNotNone(out.atoms.cell)
       self.assertFalse(out.atoms.cell is cell)
       self.assertEqual(cell.cellpar(), out.atoms.cell.cellpar())
