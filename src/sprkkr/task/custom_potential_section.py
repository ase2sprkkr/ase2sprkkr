from ..common.grammar_types import string
from ..common.options import Option

class CustomPotentialSectionDefinition(BaseValueDefintion):
  """ There is no grammar in custom potential section -
  custom sections are readed by Potential class
  """
  prefix = ''
  name_value_delimiter = '\n'

class CustomPotentialSection(Option):

      def __init__(self, name):
          super().__init__(CustomPotentialSectionDefinition(name, string))
