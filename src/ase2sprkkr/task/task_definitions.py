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
from ..common.grammar_types import mixed, flag, DefKeyword
from ..common.grammar import generate_grammar, delimitedList
from ..common.misc import lazy_value, cache
from .tasks import Task
from .outputs.kkrscf_process_output_reader import KkrScfProcessOutputReader


with generate_grammar():
  section_line_ends = pp.ZeroOrMore(pp.ZeroOrMore(pp.LineEnd().setWhitespaceChars('')) + pp.White(' \t'))

class ValueDefinition(BaseValueDefinition):
  """ This class describes the format of one value of
  a task configuration """
  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
    return pp.Suppress("=").setName('=')

  prefix = "\t"
  name_value_delimiter = '='

class SectionDefinition(BaseSectionDefinition):
  """ This class describes the format of one
  value of a task section """

  """ standard child class """
  child_class = ValueDefinition
  """ This class is used for user-added values. """
  custom_class = staticmethod(CustomOption.factory(ValueDefinition, mixed))

  """ options are delimited by newline in ouptut. """
  delimiter = '\n'
  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
      out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
      return out

  do_not_skip_whitespaces_before_name = True

class TaskDefinition(ConfDefinition):
  """ This class describes the format of a task file. """

  """ standard child class """
  child_class = SectionDefinition
  result_class = Task

  delimiter = "\n"
  @staticmethod
  @lazy_value
  def grammar_of_delimiter():
      def ws(x):
          return x.setWhitespaceChars('')
      out = (pp.Optional(section_line_ends) + pp.OneOrMore(ws(pp.LineEnd())) + pp.FollowedBy(ws(pp.Regex(r'[^\s]'))) ).suppress()
      out.setName('<newline><printable>')
      return out

  custom_class = staticmethod(CustomSection.factory(SectionDefinition))

  @classmethod
  @cache
  def custom_value_grammar(cls):
      value  = cls.child_class.custom_member_grammar()
      delim = cls.child_class.grammar_of_delimiter()
      return delimitedList(value, delim).setParseAction(lambda x: unique_dict(x.asList()))

  def __init__(self, name, sections=None,
               command='kkrscf', mpi=False, result_reader=KkrScfProcessOutputReader,
               **kwargs):
      """
      Parameters
      ---------
      ....
      command: str
        Executable to run

      mpi: bool
        Whether to run MPI version of the command

      process_class: common.process_output_reader.BaseProcessOutputReader
        Class, that runs the process and read the results
      """
      self.command = command
      self.mpi = mpi
      self.result_reader = result_reader

      super().__init__(name, sections, **kwargs)
      if not 'TASK' in self:
         self['TASK'] = SectionDefinition('TASK', [ ValueDefinition('TASK', DefKeyword(self.name)) ] )
      elif not 'TASK' in self['TASK']:
         self['TASK']['TASK'] = ValueDefinition(DefKeyword(self.name))