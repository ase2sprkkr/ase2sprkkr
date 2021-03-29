"""
This file contains the classes for definitions of Potentials:
the list of sections and their allowed (or standard) options
and their value formats. Tasks and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
"""

import functools
import pyparsing as pp
from ..common.grammar import line_end, separator as separator_grammar
from ..common.grammar_types import separator
from ..common.misc import add_to_signature
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition, \
    ConfDefinition
from .custom_potential_section import CustomPotentialSection, CustomPotentialSectionDefinition
from .potentials import Potential
from .potential_sections import PotentialSection, ASEArraySection
from ..common.misc import lazy_value

class PotValueDefinition(BaseValueDefinition):

  @staticmethod
  @lazy_value
  def _grammar_of_delimiter():
    return pp.Empty().setName(' ')

  prefix = ''
  name_value_delimiter = '\t'

  def __init__(self, *args, required=None, **kwargs):
      super().__init__(*args, required=required, **kwargs)


class Separator(PotValueDefinition):

  _counter = 0
  def __init__(self, name = None):
      if not name:
         Separator._counter += 1
         name = f'_Separator_{Separator._counter}'
      super().__init__(name, separator, is_hidden=True, name_in_grammar=False)

class PotSectionDefinition(BaseSectionDefinition):
  """ This class describes the format of one
  value of a standard potential section """

  """ The order of items in potential file is fixed """
  force_order = True
  value_name_format = '<12'

  """ standard child class """
  value_class = PotValueDefinition
  """ Adding a custom values is not allowed """
  #custom_class = staticmethod(CustomOption.factory(PotValueDefinition))
  """ A missing value is not allowed in a potential file """
  optional_value = None

  """ options are delimited by newline in ouptut. """
  delimiter = '\n'
  @staticmethod
  def _grammar_of_delimiter():
    return line_end

  def depends_on(self):
      """ The order of processing of sections during reading can be different than the order during a write. So, if the function should not be processed before given named sections, name then.

      Return
      ------
      prerequisites: [ str, str, ... ]
      """
      return []

  result_class = PotentialSection

class ASEArraySectionDefinition(PotSectionDefinition):

  @add_to_signature(PotSectionDefinition.__init__)
  def __init__(self, *args, array_name, **kwargs):
      super().__init__(*args, **kwargs)
      self.array_name = array_name

  def depends_on(self):
      """ Array size is required """
      return [ 'SITES' ]

  result_class = ASEArraySection

class PotentialDefinition(ConfDefinition):
  """ This class describes the format of a potential file """

  """ standard child class """
  section_class = PotSectionDefinition
  result_class = Potential


  """ The order of items in potential file is fixed """
  force_order = True

  delimiter="*"*80 + "\n"

  @staticmethod
  @lazy_value
  def _grammar_of_delimiter():
      return pp.And([
                     pp.OneOrMore(line_end),
                     separator_grammar,
                     pp.OneOrMore(line_end)
             ]).setName('*'*80 + '<newline>')

  custom_class = CustomPotentialSection

  @staticmethod
  @lazy_value
  def _custom_section_value():
      return pp.SkipTo(separator_grammar | pp.StringEnd()).setName('SkipTo ******|EndOfFile')
