""" Custom potential sections are the parts of Potential, that whose content is not parsed.

I.e. these sections can has any content (they are readed up to the section separator).
"""


from ..common.grammar_types import BaseType
from ..common.options import CustomOption
from ..common.configuration_definitions import BaseValueDefinition
from ..common.misc import class_property, cache
import pyparsing as pp
from ..common.grammar import separator, line_end
from ..common.misc import lazy_value
import re
import pyparsing

class CustomPotentialSectionDefinition(BaseValueDefinition):
  """ There is no grammar in a custom potential section -
  custom sections are readed by Potential class
  """
  mandatory = False
  """ Obviously, the custom sections are not required """
  prefix = ''
  name_value_delimiter = '\n'
  """ The content of the section is delimited from the name by a newline """

class CustomSectionToken(pp.Token):
   """ The grammar for a custom section - i.e. for unknown section, whose
   content is let as is.

   The grammar just reads all up to the section separator.
   """

   pattern = re.compile('\n' + separator.pattern + '[ \r\t]*\n',  re.DOTALL)
   name = 'CustomSection'

   def parseImpl(self, instr, loc, doActions = True):
       result = self.pattern.search(instr,loc)
       if result:
          out = instr[loc:result.start()]
          loc = result.start()
       else:
          out = instr[loc:]
          loc = len(instr)
       return loc, pp.ParseResults(out.strip())

class SectionString(BaseType):
      """
      The grammar_type of a custom section - i.e. string, that
      ends with a section separator.

      This grammar_type as used as a value type for the custom section.
      """

      delimiter_pattern = '(?:[ \t\r]*(?:\n[ \t\r]*)*)*\n' +separator.pattern + '(?:[ \t\r]*(?:\n[ \t\r]*))*\n'

      @staticmethod
      @lazy_value
      def grammar_of_delimiter():
          return pp.Regex(SectionString.delimiter_pattern).setName('*'*79 + '<newline>').suppress()

      @class_property
      @cache
      def _grammar(cls):
          return CustomSectionToken()

      def grammar_name(self):
          return '<all up to the end of the section>'

      def write(self, f, value):
          super().write(f, value)
          f.write('\n')

SectionString.I = SectionString()

class CustomPotentialSection(CustomOption):
      """
      Unknown sections of the potential file are mapped to a "section"
      of this type.

      In fact, it is not a Section - a container - but just an Option,
      that holds a string value: a content of the section.
      """
      def __init__(self, name, container=None):
          super().__init__(CustomPotentialSectionDefinition(name, SectionString.I), container)

      def _depends_on(self):
          return []

      def _set_from_atoms(self, atoms, io_data):
          pass

      def _update_atoms(self, atoms, io_data):
          pass

      def reset(self):
          return self.remove()
