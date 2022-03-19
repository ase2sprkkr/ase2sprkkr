if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..definitions.sections.lattice import LatticeSection, LatticeSectionDefinition
from ..definitions.sections.sites import SitesSectionDefinition
from ..potentials import Potential
from ..potential_definitions import PotentialDefinition
from ase.build import bulk
import io

class TestStructure(TestCase):

  def test_lattice(self):
    lsd=LatticeSectionDefinition('LATTICE')
    ssd=SitesSectionDefinition('SITES')
    pd=PotentialDefinition('pot', [lsd, ssd])
    pot=Potential(definition = pd)
    sec=pot.LATTICE

    for e in ['Co', 'Ni', 'U', 'Ge', 'As', 'Na', 'In', 'Pr']:
       atoms = bulk(e)
       pot.atoms = atoms
       cell = atoms.cell
       s = io.StringIO()
       pot.save_to_file(s)
       s.seek(0)
       out = pd.read_from_file(s)

       a=out.LATTICE.as_dict()
       b=sec.as_dict()

       self.assertEqual(out.LATTICE.as_dict(), sec.as_dict())
       self.assertIsNotNone(out.atoms)
       new_cell = out.atoms.cell
       self.assertIsNotNone(new_cell)
       self.assertFalse(cell is new_cell)
       self.assertEqual(cell * 1, new_cell * 1)
       self.assertEqual(out.atoms.positions, atoms.positions)
