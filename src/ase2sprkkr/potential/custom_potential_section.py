from ..common.grammar_types import BaseType
from ..common.options import Option
from ..common.configuration_definitions import BaseValueDefinition
from ..common.misc import class_property, cache
import pyparsing as pp
from ..common.grammar import separator, line_end
from ..common.misc import lazy_value
import re
import pyparsing

class CustomPotentialSectionDefinition(BaseValueDefinition):
  """ There is no grammar in custom potential section -
  custom sections are readed by Potential class
  """
  prefix = ''
  name_value_delimiter = '\n'

class CustomSectionToken(pp.Token):
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

      delimiter_pattern = '(?:[ \t\r]*(?:\n[ \t\r]*)*)*\n' +separator.pattern + '(?:[ \t\r]*(?:\n[ \t\r]*))*\n'

      @staticmethod
      @lazy_value
      def grammar_of_delimiter():
          return pp.Regex(SectionString.delimiter_pattern).setName('*'*80 + '<newline>').suppress()

      @class_property
      @cache
      def _grammar(cls):
          return CustomSectionToken()

      def grammar_name(self):
          return '<all up to end of section>'

SectionString.I = SectionString()

class CustomPotentialSection(Option):
      def __init__(self, name, container=None):
          super().__init__(CustomPotentialSectionDefinition(name, SectionString.I), container)

      def _depends_on(self):
          return []

      def _set_from_atoms(self, atoms, io_data):
          pass

      def _update_atoms(self, atoms, io_data):
          pass