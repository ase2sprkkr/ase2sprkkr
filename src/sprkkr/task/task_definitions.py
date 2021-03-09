"""
This file contains the classes for definitions of Tasks:
the list of sections and their allowed (or standard) options
and their value formats. Tasks and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
"""

import functools
import pyparsing as pp
from ..common.grammar import optional_line_end, line_end
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition, \
    ConfDefinition, \
    unique_dict
from ..common.options import CustomOption
from ..common.conf_containers import CustomSection
from ..common.grammar_types import mixed, flag
from ..common.misc import lazy_value
from .tasks import Task

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
  def _grammar_of_delimiter():
   return optional_line_end


class TaskDefinition(ConfDefinition):
  """ This class describes the format of a task file. """

  """ standard child class """
  section_class = SectionDefinition
  result_class = Task

  delimiter = "\n"
  @lazy_value
  @staticmethod
  def _grammar_of_delimiter():
    return line_end + pp.OneOrMore(line_end)

  custom_class = staticmethod(CustomSection.factory(SectionDefinition))
  @lazy_value
  @staticmethod
  def _custom_section_value():
      value  = SectionDefinition._grammar_of_delimiter() + SectionDefinition.custom_value()
      return pp.OneOrMore(value).setParseAction(lambda x: unique_dict(x.asList()))
