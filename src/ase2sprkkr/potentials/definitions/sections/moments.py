from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection as PotSection, AtomicTypePotentialSection
from ....common.grammar_types import NumpyArray
from ....common.configuration_definitions import SeparatorDefinition
import re


class MomentSection(PotSection):
    pass


class MomentsSection(AtomicTypePotentialSection):

    property_name = 'moments'

    def write_data(self, typ, section, index):
        section.TYPE = index
        section.DATA = typ.moments.as_tuple()


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
                       name_regex=re.compile(r'MOMENTS(?:(?: *[A-Z]{3}){5})'),
                       write_alternative_name=True,
                       )

  result_class = MomentSection
  repeated_class = MomentsSection


section = MomentsSectionDefinition
