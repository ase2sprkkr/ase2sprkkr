"""
This file contains the classes for definitions of parts of output
files, from that the result data are obtained.

Although the files are used only for reading, the same
grammar definition approach as for potential and input parameters
files is used
"""

import functools
import pyparsing as pp
from ..common.grammar import line_end
from ..common.configuration_definitions import \
    BaseValueDefinition, \
    BaseSectionDefinition
from ..common.misc import lazy_value

class OutputValueDefinition(BaseValueDefinition):
  """ Value in an output file, of a form 'NAME   VALUE' """

  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
    return pp.WordStart()

  prefix = ''

  def __init__(self, *args, required=True, **kwargs):
      super().__init__(*args, required=required, **kwargs)

class OutputValueEqualDefinition(OutputValueDefinition):
  """ Value in an output file, of a form 'NAME=VALUE' (spaces possible) """
  name_value_delimiter = '='

  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
    return pp.Suppress("=").setName('=')

class OutputNonameValueDefinition(OutputValueDefinition):
  """ Value in an output file, that has no name, there is just the value
  (identified by its position after some other known stuff)
  """

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
