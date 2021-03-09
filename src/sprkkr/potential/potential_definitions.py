"""
This file contains the classes for definitions of Potentials:
the list of sections and their allowed (or standard) options
and their value formats. Tasks and potentials have a simlilar 
structure, so they share common functionalities from 
sprkkr.common.configuration_definitions
"""

import functools
import pyparsing as pp
from grammar import line_end, separator
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition, \
    ConfDefinition
from .import CustomPotentialSection, CustomPotentialSectionDefinition
from ..common.misc import lazy_value

class PotValueDefinition(BaseValueDefinition):

  @lazy_value
  @staticmethod
  def _grammar_of_delimiter():
    return pp.Empty().setName(' ')

  prefix = ''
  name_value_delimiter = '\t'


class PotDataDefinition(PotValueDefinition):
  name_in_grammar = False

class Separator(PotValueDefinition):
  
  name_in_grammar = False
 
  _counter = 0
  def __init__(self, name = None):
      if not name:
         Separator.counter += 1
         name = '_separator_{Separator.counter}'
      super.__init__(name, separator)

class PotSectionDefinition(BaseSectionDefinition):
  """ This class describes the format of one 
  value of a standard potential section """

  """ The order of items in potential file is fixed """
  force_order = True

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


class PotentialDefinition(ConfDefinition):
  """ This class describes the format of a potential file """

  """ standard child class """
  section_class = PotSectionDefinition
  result_class = Potential 


  """ The order of items in potential file is fixed """
  force_order = True

  @lazy_value
  @staticmethod
  def _grammar_of_delimiter():
      return (pp.OneOrMore(line_end) + separator + pp.OneOrMore(line_end)).setName('*'*80 + '\n')

  custom_class = CustomPotentialSection
  @lazy_value
  @staticmethod
  def _custom_section_value(section_names):
      return pp.SkipTo(separator | pp.StringEnd())
