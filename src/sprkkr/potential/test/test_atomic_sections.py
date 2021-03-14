from .init_tests import TestCase, patch_package
__package__ = patch_package(__package__)

import os
from ..lattice_section import LatticeSection, LatticeSectionDefinition
from ..sites_section import SitesSectionDefinition
from ..potentials import Potential
from ..potential_definitions import PotentialDefinition
from ase.build import bulk
import io

def add_section_create_potential(definition, pot, name, cls):
    definition[name] = cls(name)
    pot._add(definition[name].create_object(pot))
    s = io.StringIO()
    pot.save_to_file(s)
    s.seek(0)
    pot2 = Potential(definition = definition)
    pot2.read_from_file(s)
    return pot2


class TestLattice(TestCase):

  def test_lattice(self):
    lsd=LatticeSectionDefinition('LATTICE')
    pd=PotentialDefinition('pot', [lsd])
    pot=Potential(definition = pd)
    sec=pot.LATTICE

    for e in ['Co', 'Ni', 'U', 'Ge', 'As', 'Na', 'In', 'Pr']:
       atoms = bulk(e)
       pot.atoms = atoms
       cell = atoms.cell
       s = io.StringIO()
       sec.save_to_file(s)
       s.seek(0)
       out = pd.read_from_file(s)

       a=out.LATTICE.to_dict()
       b=sec.to_dict()

       self.assertEqual(out.LATTICE.to_dict(), sec.to_dict())
       self.assertIsNone(out.atoms)
       out._atoms_io_data.set_atoms_positions([])
       self.assertIsNotNone(out.atoms)
       new_cell = out.atoms.cell
       self.assertIsNotNone(new_cell)
       self.assertFalse(cell is new_cell)
       self.assertEqual(cell.cellpar(), new_cell.cellpar())
       #The resulting cell could be different!!!
       #self.assertEqual(cell*1., new_cell*1.)

    other = bulk(e)

    pot2 = add_section_create_potential(pd, pot, 'SITES', SitesSectionDefinition)
    self.assertFalse(pot.atoms is pot2.atoms)
    self.assertEqual(pot.atoms.positions, pot2.atoms.positions)
