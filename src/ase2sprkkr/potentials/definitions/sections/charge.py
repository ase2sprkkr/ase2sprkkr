from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection as PotSection, RepeatedPotentialSection
from ....common.grammar_types import NumpyArray, RawData
from ....common.unique_values import UniqueValuesMapping
from ....common.warnings import DataValidityWarning
from ....common.configuration_definitions import SeparatorDefinition
import re


class ChargeSection(PotSection):
  pass


class ChargesSection(RepeatedPotentialSection):

  def _set_from_atoms(self, atoms, write_io_data):
      self.clear()
      if not write_io_data.has_converged_data(self._container):
          return
      for site, i in write_io_data.sites.unique_items():
          charge = self.add(i)
          charge.TYPE = i
          charge.DATA = site.charge.raw_value

  def _update_atoms(self, atoms, read_io_data):
      if len(self):
          try:
              for site, id in UniqueValuesMapping.from_values(atoms.sites).unique_items():
                  site.charge = self[id - 1]['DATA']()
          except KeyError:
              DataValidityWarning("There is no CHARGE data for the atom type with id {id}.")

  def _depends_on(self):
      return 'TYPES'


class ChargeSectionDefinition(PotSectionDefinition):

  def __init__(self, name='CHARGE', **kwargs):
      V = PotValueDefinition
      members = [
          V('TYPE', int),
          V('DATA', NumpyArray(line_length=100, shape=(2,-1), ends_with=re.compile("\n(={79}|-{79})"),
                               ends_with_str = "=" * 79, item_format='% .14E', indented=1),
                    name_in_grammar=False,
           ),
          V('FULLPOT', RawData(ends_with=re.compile("\n?={79}"), ends_with_str="=" * 79, condition=lambda x: len(x)), required=False,
            name_in_grammar=False, write_condition = lambda x: x() != ''),
          SeparatorDefinition('=', length=79)
      ]
      super().__init__(name, members, has_hidden_members=True, is_repeated=True, is_optional=True)

  result_class = ChargeSection
  repeated_class = ChargesSection


section = ChargeSectionDefinition
