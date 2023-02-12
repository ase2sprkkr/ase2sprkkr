from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import UniqueListSection

from ....common.grammar_types import Table
from ....sprkkr.reference_systems import ReferenceSystem

class ReferenceSystemSection(UniqueListSection):
  _value_name = 'reference_systems'
  _value_class = ReferenceSystem

  def _set_from_atoms(self, atoms, io_data):
      super()._set_from_atoms(atoms, io_data)
      self['NREF'].set(io_data.reference_systems.len_of_unique())

class ReferenceSystemSectionDefinition(PotSectionDefinition):

  def __init__(self, name='REFERENCE SYSTEM',
                     alternative_names='REFERENCE SYSTEM FOR TIGHT BINDING MODE',
                     write_alternative_name=True,
                     **kwargs):
      V = PotValueDefinition
      members = [
          V('NREF', int),
          V('DATA', Table({'VREF': float, 'RMTREF' :float}, numbering='IREF')),
      ]
      super().__init__(name, members, has_hidden_members=True, alternative_names = alternative_names)

  def validate(self, values, why='set'):
      return values['NREF'] == len(values['DATA']) or \
             'Number of reference systems differs from NREF'

  result_class = ReferenceSystemSection

section = ReferenceSystemSectionDefinition
