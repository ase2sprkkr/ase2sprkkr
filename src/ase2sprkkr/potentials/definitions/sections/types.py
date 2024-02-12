from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import UniqueListSection
from ....common.grammar_types import Table
from ....sprkkr.atomic_types import AtomicType


class TypesSection(UniqueListSection):
  _value_name = 'types'
  _value_class = AtomicType


class TypesSectionDefinition(PotSectionDefinition):

  def __init__(self, name='TYPES', **kwargs):
      V = PotValueDefinition
      members = [
          V('DATA', Table({'TXT': str, 'ZT' :int, 'NCORT': int, 'NVALT': int, 'NSEMCORSHLT': int},
                          numbering='IT', free_header=True)
          ),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = TypesSection


section = TypesSectionDefinition
