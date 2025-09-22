from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection as PotSection, AtomicTypePotentialSection
from ....common.grammar_types import NumpyArray, RawData
from ....common.configuration_definitions import SeparatorDefinition
import re


class PotentialSection(PotSection):
    pass


class PotentialsSection(AtomicTypePotentialSection):

    property_name = 'potential'
    property_label = 'radial potential'


class PotentialSectionDefinition(PotSectionDefinition):

  def __init__(self, name='POTENTIAL', **kwargs):
      V = PotValueDefinition
      members = [
          V('TYPE', int),
          V('DATA', NumpyArray(line_length=100, shape=(2,-1), ends_with=re.compile("\n?(={79}|-{79}|NFP)"),
                               item_format='% .14E', indented=1),
                    name_in_grammar=False,
           ),

          V('FULLPOT', RawData(ends_with=re.compile("\n?(={79}|-{79})"),
                               condition = lambda x: x != ''
                               ),
                    is_required=False,
                    name_in_grammar=False,
                    write_condition = lambda x: x() != ''
          ),
          SeparatorDefinition('=', length=79)
      ]
      super().__init__(name, members, has_hidden_members=True, is_repeated=True, is_optional=True)

  result_class = PotentialSection
  repeated_class = PotentialsSection


section = PotentialSectionDefinition
