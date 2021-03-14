from .potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from .atoms_sections import AtomsSection

from  ..common.grammar_types import DefKeyword, Array, Table, Integer

class SitesSection(AtomsSection):
  """ This section retrieves the atomic positions and
      it creates (during reading) the ASE Atoms object """

  def _set_from_atoms(self):
      aoid=self._atoms_io_data
      self['SCALED_ATOMIC_POSITIONS'].set(
        aoid.bravais_cell.scaled_positions(
          aoid.atoms.positions
        )
      )

  def _process(self):
      sap = self['SCALED_ATOMIC_POSITIONS']()
      if sap is not None:
          aiod = self._atoms_io_data
          pos = aiod.bravais_cell.cartesian_positions(sap)
          aiod.set_atoms_positions(pos)

class SitesSectionDefinition(PotSectionDefinition):

  def __init__(self, name='SITES', **kwargs):
      V = PotValueDefinition
      members = [
          V('CARTESIAN', bool, fixed_value=True),
          V('BASSSCALE', Array(float, length=3), fixed_value=[1.,1.,1.]),
          V('SCALED_ATOMIC_POSITIONS', Table({'QX': float, 'QY' : float, 'QZ': float}, numbering='IQ')),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = SitesSection
