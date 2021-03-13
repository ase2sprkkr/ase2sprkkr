from .potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from .atom_sections import AtomSection

from  ..common.grammar_types import DefKeyword, Sequence, Table, Integer

class SitesSection(AtomSection):

  def _set_from_atoms(self):
      self['SCALED_ATOMIC_POSITIONS'].set(
        self.atoms.get_scaled_positions()
      )

  def _process(self):
      sap = self['SCALED_ATOMIC_POSITIONS']()
      if sap is not None:
          aiod = self._atoms_io_data
          pos = aiod.cell.cartesian_positions(sap)
          aiod.set_new_positions(pos)

class SitesSectionDefinition(PotSectionDefinition):

  def __init__(self, name, **kwargs):
      V = PotValueDefinition
      members = [
          V('CARTESIAN', bool, fixed_value=True),
          V('BASSSCALE', Array(float, length=3), fixed_value=[1.,1.,1.]),
          V('SCALED_ATOMIC_POSITIONS', Table({'QX': float, 'QY' : float, 'QZ': float}, numbering='IQ')),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = SitesSection
