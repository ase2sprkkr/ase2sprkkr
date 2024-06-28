if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..definitions.sections.lattice import LatticeSectionDefinition  # NOQA: E402
from ..definitions.sections.sites import SitesSectionDefinition      # NOQA: E402
from ..potentials import Potential                                   # NOQA: E402
from ..potential_definitions import PotentialDefinition              # NOQA: E402
from ase.build import bulk                                           # NOQA: E402
import io                                                            # NOQA: E402


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

       with self.almost_equal_precision(decimal=12):
           self.assertEqual(out.LATTICE.as_dict(), sec.as_dict())
       assert out.atoms is not None
       new_cell = out.atoms.cell
       assert new_cell is not None
       self.assertFalse(cell is new_cell)
       self.assertEqual(cell * 1, new_cell * 1)
       self.assertEqual(out.atoms.positions, atoms.positions)
