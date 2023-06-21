from ..sprkkr.configuration import ConfigurationValueDefinition, ConfigurationFileDefinition, \
                                   CustomConfigurationValue

from ..common.decorators import cached_class_property, cache
from ..common.grammar_types  import unsigned, Array, Table, RestOfTheFile, Keyword, GrammarType, \
                                     pot_mixed, line_string, line_end, Separator as GTSeparator
from ..common.grammar import generate_grammar
import pyparsing as pp
from ..common.decorators import cached_class_property
import sys
from .output_files import OutputFile

class OutputFileValueDefinition(ConfigurationValueDefinition):
  """ This class describes the format of one value of
  a header of an output file """
  @cached_class_property
  def grammar_of_delimiter():
    return pp.Empty().setName(' ')

  prefix = " "
  name_value_delimiter = '\t'

  type_from_type_map = { str: line_string }
  type_of_dangerous = pot_mixed

class Separator(OutputFileValueDefinition):
  """
  A special (hidden) value, that appears in a output file header

  The separator is a blank line
  """
  _counter = 0
  def __init__(self, name = None):
      if not name:
         Separator._counter += 1
         name = f'_Separator_{Separator._counter}'
      with generate_grammar():
        super().__init__(name, GTSeparator(pp.Empty()), is_hidden=True, name_in_grammar=False)


class OutputFileDefinition(ConfigurationFileDefinition):
  """ This class describes the format of one
  value of a standard potential section """

  force_order = True
  """ The order of items in potential file is fixed """

  value_name_format = '<12'

  child_class = OutputFileValueDefinition
  """ standard child class """

  custom_class = None
  """ No custom members in the output files """

  delimiter = '\n'
  """ options are delimited by newline in ouptut. """

  @staticmethod
  def grammar_of_delimiter():
    return line_end

  result_class = OutputFile

  configuration_type_name = 'OUTPUT FILE'

@cache
def output_file_header():
    """ Return the members of the common output file header, up to
    the first ``KEYWORD`` value (which differs for each output file)"""
    V = OutputFileValueDefinition
    return [
      V('TITLE', str),
      V('SYSTEM', str),
      V('NQ_EFF', int),
      V('NT_EFF', int),
      Separator(),
      V('NE', int),
      V('IREL', int),
      V('EFERMI', float),
      V('INFO', str),
      Separator(),
      V('ORBITALS', Table({'NLQ' : unsigned}, numbering='IQ', flatten=True), name_in_grammar = False),
      V('TYPES', Table({'TXT_T': str, 'CONC': float, 'NAT': int, 'IQAT': Array(int)}, numbering='IT'), name_in_grammar=False),
    ]

def create_output_file_definition(keyword, add, name=None,
                                  cls=OutputFileDefinition,
                                  **kwargs):
    """ Create a definition of a output file. One should supply the additional values
    containing the actual data.
    """
    V = OutputFileValueDefinition
    if isinstance(keyword, str):
        if name is None:
           name = keyword
        keyword = V('KEYWORD', Keyword(keyword) )
    elif isinstance(keyword, GrammarType):
        keyword = V('KEYWORD', keyword )
    fields = [ keyword ]
    fields.extend( output_file_header() )
    if add:
       fields.extend(add)
    return cls(f'{name}', fields, **kwargs)
