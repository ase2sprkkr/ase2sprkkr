from ..potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ..potential_sections import PotentialSection

from  ...common.grammar_types import DefKeyword, Array, Table, Integer
from ...ase.sprkkr_atoms import SprKkrAtoms

class SitesSection(PotentialSection):
  """ This section retrieves the atomic positions and
      it creates (during reading) the ASE Atoms object """

  def _set_from_atoms(self, atoms, write_io_data):
      self['SCALED_ATOMIC_POSITIONS'].set(
          atoms.positions / ( write_io_data['lattice.alat'] * self['BASSCALE']())
      )

  def _update_atoms(self, atoms, read_io_data):
      positions = self['SCALED_ATOMIC_POSITIONS']() * \
            (read_io_data['lattice.alat'] * self['BASSCALE']())
      new = True
      try:
         if atoms:
            atoms.set_positions(positions)
            new = False
      except ValueError:
            pass
      if new:
         atoms = SprKkrAtoms(positions = positions, cell = read_io_data['lattice.cell'])
         return atoms

class SitesSectionDefinition(PotSectionDefinition):

  def __init__(self, name='SITES', **kwargs):
      V = PotValueDefinition
      members = [
          V('CARTESIAN', bool, fixed_value=True),
          V('BASSCALE', Array(float, length=3), fixed_value=[1.,1.,1.]),
          V('SCALED_ATOMIC_POSITIONS', Table({'QBAS(X)': float, 'QBAS(Y)' : float, 'QBAS(Z)': float}, numbering='IQ',free_header=True)),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = SitesSection
