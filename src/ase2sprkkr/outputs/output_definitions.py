"""
This file contains the classes for definitions of parts of output
files. Although the files are used only for reading, the same
grammar as for potential and input files is used
"""

import functools
import pyparsing as pp
from ..common.grammar import line_end
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition
from ..common.misc import lazy_value

class OutputValueDefinition(BaseValueDefinition):

  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
    return pp.WordStart()

  prefix = ''
  name_value_delimiter = '='

  def __init__(self, *args, required=True, **kwargs):
      super().__init__(*args, required=required, **kwargs)

class OutputValueEqualDefinition(OutputValueDefinition):

  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
    return pp.Suppress("=").setName('=')

class OutputNonameValueDefinition(OutputValueDefinition):

  name_in_grammar = False

  def __init__(self, *args, required=True, **kwargs):
      super().__init__(*args, required=required, **kwargs)



class OutputSectionDefinition(BaseSectionDefinition):
  """ This class describes the format of one
  value of a standard potential section """

  """ The order of items in output file is fixed """
  force_order = True
  name_in_grammar = False

  """ standard child class """
  child_class = OutputValueDefinition
  custom_class = None

  """ options are delimited by newline in ouptut. """
  delimiter = '\n'
  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
      out = (pp.Optional(line_end) + pp.WordStart()).suppress()
      return out
