from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import ASEArraySection

from ....common.grammar_types import Table


class ChargeMomentsSection(ASEArraySection):

    pass


class ChargeMomentsSectionDefinition(PotSectionDefinition):

  array_name = 'charge_moments'

  def __init__(self, name='CHARGE MOMENTS', **kwargs):
      V = PotValueDefinition
      members = [
          V('DATA', Table(numbering='IQ', grouping=True, group_size='NLMTOP-CHR',
                          CMNTQ=float, numbering_format='{:>10}', grouping_format='{:>3}',
                          flatten=True),
                    name_in_grammar=False, is_optional=True
           )
      ]
      super().__init__(name, members, has_hidden_members=True, is_optional=True)

  result_class = ChargeMomentsSection


section = ChargeMomentsSectionDefinition