"""
This file contains the classes for definitions of InputParameterss:
the list of sections and their allowed (or standard) options
and their value formats. InputParameterss and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
"""

import pyparsing as pp
from ..common.container_definitions import dict_from_parsed
from ..sprkkr.configuration import \
         ConfigurationValueDefinition, ConfigurationSectionDefinition, ConfigurationFileDefinition, \
         CustomConfigurationValue, CustomConfigurationSection
from ..common.grammar_types import mixed, flag, DefKeyword
from ..common.grammar import generate_grammar, delimitedList
from ..common.decorators import cached_class_property, cache
from .input_parameters import InputParameters, InputSection

with generate_grammar():
  section_line_ends = pp.ZeroOrMore(pp.ZeroOrMore(pp.LineEnd().setWhitespaceChars('')) + pp.White(' \t'))


class InputValueDefinition(ConfigurationValueDefinition):
  """ This class describes the format of one value of
  a task configuration """
  @cached_class_property
  def grammar_of_delimiter():
    return pp.Suppress("=").setName('=')

  prefix = "\t"
  name_value_delimiter = '='

  type_from_type_map = { bool : flag }
  type_of_dangerous = mixed


class InputSectionDefinition(ConfigurationSectionDefinition):
  """ This class describes the format of one
  value of a task section """

  child_class = InputValueDefinition
  """ standard child class """

  result_class = InputSection
  """ The standard class for InputParameters section """

  custom_class = staticmethod(CustomConfigurationValue.factory(InputValueDefinition, mixed))
  """ Factory for custom values in the input sections. """

  delimiter = '\n'
  """ options are delimited by newline in ouptut. """

  @cached_class_property
  def grammar_of_delimiter():
      out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
      return out

  do_not_skip_whitespaces_before_name = True


class InputParametersDefinition(ConfigurationFileDefinition):
  """ This class describes the format of a task file. """

  save_hook = None
  """ Input parameters can have a save_hook defined, which is executed before saving the
  parameters. The arguments of the hook are ``filename, atoms, self``. """

  child_class = InputSectionDefinition
  """ Sections of the :class:`InputParameters` are defined by :class:`InputSectionDefinition` """

  result_class = InputParameters
  """ The parsing of a potential file results in an instance of :class:`InputParameters` """

  custom_class = staticmethod(CustomConfigurationSection.factory(InputSectionDefinition))
  """ The class factory for custom sections in the container """

  configuration_type_name = 'INPUT PARAMETERS'
  """ Name of the container type in the runtime documentation """

  delimiter = "\n"
  """ Sections are delimited by newline in the output """

  @cached_class_property
  def grammar_of_delimiter():
      def ws(x):
          return x.setWhitespaceChars('')
      out = (pp.Optional(section_line_ends) + pp.OneOrMore(ws(pp.LineEnd())) + pp.FollowedBy(ws(pp.Regex(r'[^\s]'))) ).suppress()
      out.setName('<newline><printable>')
      return out

  @classmethod
  @cache
  def custom_value_grammar(cls):
      value = cls.child_class.custom_member_grammar()
      delim = cls.child_class.grammar_of_delimiter()
      return delimitedList(value, delim).\
            setParseAction(lambda x: dict_from_parsed(x.asList(), lambda x: True))

  def _generic_info(self):
      return f"Input parameters for task {self.name}"

  def __init__(self, name, members=None,
               executable='kkrscf', mpi=True, result_reader=None,
               **kwargs):
      """
      Parameters
      ---------
      others:
        For the meaning of the others parameters please see :class:`ConfigurationRootDefinition`

      executable: str
        Executable to run

      mpi: bool
        Whether to run MPI version of the executable

      result_reader: common.process_output_reader.ProcessOutputReader
        Class, that runs the process and read the results. Default NONE
        means, that the class is determined from the TASK name
        (see InputParameters.result_reader)
      """
      self.executable = executable
      self.mpi = mpi
      self.result_reader = result_reader

      super().__init__(name, members, **kwargs)
      if not 'TASK' in self:
         self['TASK'] = InputSectionDefinition('TASK', [ InputValueDefinition('TASK', DefKeyword(self.name), name_in_grammar=False) ] )
      elif not 'TASK' in self['TASK']:
         self['TASK']['TASK'] = InputValueDefinition(DefKeyword(self.name), name_in_grammar=False)
