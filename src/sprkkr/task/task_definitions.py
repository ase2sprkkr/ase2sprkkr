"""
This file contains the classes for definitions of Tasks:
the list of sections and their allowed (or standard) options
and their value formats. Tasks and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
"""

import functools
import pyparsing as pp
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition, \
    ConfDefinition, \
    unique_dict
from ..common.options import CustomOption
from ..common.conf_containers import CustomSection
from ..common.grammar_types import mixed, flag, Keyword
from ..common.grammar import generate_grammar
from ..common.misc import lazy_value
from .tasks import Task

with generate_grammar():
  section_line_ends = pp.ZeroOrMore(pp.ZeroOrMore(pp.LineEnd().setWhitespaceChars('')) + pp.White(' \t'))

class ValueDefinition(BaseValueDefinition):
  """ This class describes the format of one value of
  a task configuration """
  @staticmethod
  @lazy_value
  def _grammar_of_delimiter():
    return pp.Suppress("=").setName('=')

  prefix = "\t"
  name_value_delimiter = '='

class SectionDefinition(BaseSectionDefinition):
  """ This class describes the format of one
  value of a task section """

  """ standard child class """
  value_class = ValueDefinition
  """ This class is used for user-added values. """
  custom_class = staticmethod(CustomOption.factory(ValueDefinition))
  """ Missing value means the flag type """
  optional_value = flag

  """ options are delimited by newline in ouptut. """
  delimiter = '\n'
  @staticmethod
  @lazy_value
  def _grammar_of_delimiter():
      out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
      return out


class TaskDefinition(ConfDefinition):
  """ This class describes the format of a task file. """

  """ standard child class """
  section_class = SectionDefinition
  result_class = Task

  delimiter = "\n"
  @staticmethod
  @lazy_value
  def _grammar_of_delimiter():
      def ws(x):
          return x.setWhitespaceChars('')
      out = (pp.Optional(section_line_ends) + pp.OneOrMore(ws(pp.LineEnd())) + pp.FollowedBy(ws(pp.Regex(r'[^\s]'))) ).suppress()
      out.setName('<newline><printable>')
      return out

  custom_class = staticmethod(CustomSection.factory(SectionDefinition))
  @lazy_value
  @staticmethod
  def _custom_section_value():
      value  = SectionDefinition._grammar_of_delimiter() + SectionDefinition.custom_value()
      return pp.OneOrMore(value).setParseAction(lambda x: unique_dict(x.asList()))

  def __init__(self, name, sections=None, **kwargs):
      super().__init__(name, sections, **kwargs)
      if not 'TASK' in self:
         self['TASK'] = SectionDefinition('TASK', [ ValueDefinition('TASK', Keyword(self.name)) ] )
      elif not 'TASK' in self['TASK']:
         self['TASK']['TASK'] = ValueDefinition(Keyword(self.name))
