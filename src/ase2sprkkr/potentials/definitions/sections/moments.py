from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection as PotSection, RepeatedPotentialSection
from ....common.grammar_types import NumpyArray
from ....common.unique_values import UniqueValuesMapping
from ....common.warnings import DataValidityWarning
from ....common.configuration_definitions import SeparatorDefinition


class MomentSection(PotSection):
  pass


class MomentsSection(RepeatedPotentialSection):

  def _set_from_atoms(self, atoms, write_io_data):
      self.clear()
      if not write_io_data.has_converged_data(self._container):
          return
      for site, i in write_io_data.sites.unique_items():
          charge = self.add(i)
          charge.TYPE = i
          charge.DATA = site.moments.as_tuple()

  def _update_atoms(self, atoms, read_io_data):
      if len(self):
          try:
              for site, id in UniqueValuesMapping.from_values(atoms.sites).unique_items():
                  site.moments = self[id - 1]['DATA']()
          except KeyError:
              DataValidityWarning("There is no CHARGE data for the atom type with id {id}.")

  def _depends_on(self):
      return 'TYPES'


class MomentsSectionDefinition(PotSectionDefinition):

  def __init__(self, name='MOMENTS', **kwargs):
      V = PotValueDefinition
      members = [
          V('TYPE', int),
          V('DATA', NumpyArray(lines=1, shape=(5,), written_shape=(1,5), item_format='% .14E', indented=1), name_in_grammar=False),
          SeparatorDefinition('=', length=79)
      ]
      super().__init__(name, members, has_hidden_members=True, is_repeated=True, is_optional=True,
                       alternative_names='MOMENTS        QEL  NOS  SMT  OMT  HFF',
                       write_alternative_name=True,
                       )

  result_class = MomentSection
  repeated_class = MomentsSection


section = MomentsSectionDefinition
