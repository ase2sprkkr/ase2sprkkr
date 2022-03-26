"""
Various pyparsing grammar elements and a few useful routines.
"""

from contextlib import contextmanager
import pyparsing as pp

@contextmanager
def generate_grammar():
    """ Set the pyparsing newline handling
        and then restore to the original state """
    try:
      old = None
      pe = pp.ParserElement
      if hasattr(pe, "DEFAULT_WHITE_CHARS"):
          old = pe.DEFAULT_WHITE_CHARS
      pe.setDefaultWhitespaceChars(' \t\r')
      yield
    finally:
      if old is not None:
        pe.setDefaultWhitespaceChars(old)

def replace_whitechars(expr):
    expr = expr.copy()
    expr.setWhitespaceChars(' \t\r')
    return expr


with generate_grammar():
  optional_line_end = pp.Suppress(pp.LineEnd() | pp.WordStart() ).setName(' ')
  """ Grammar for an optinal newline """
  line_end = pp.Suppress(pp.LineEnd()).setName('\n')
  """ Grammar for a required newline """
  end_of_file = (pp.Regex(r'[\s]*') + pp.StringEnd()).suppress().setName('<EOF>')
  """ Grammar for an end of file (ending whitespaces are allowed) """

  separator_pattern = r'\*'*10+r'\**'
  """ Pattern for separating sections in an input file """
  separator = pp.Regex(separator_pattern).setName("**********[***....]").suppress()
  """ Grammar for separating sections in an input file """
  separator.pattern = separator_pattern

  optional_quote = pp.Optional("'").suppress()
  """ Grammar for an optional quote """

def delimitedList(expr, delim):
  """ Delimited list with already suppressed delimiter (or with a in-results-wanted one) """
  return expr + pp.ZeroOrMore(delim + expr)

def addConditionEx(self, condition, message):
  """ Add check condition to the pyparsing ParseElement,
  that, if it failed, raise a parse exception with a given message. """

  def check_condition(s, loc, tocs):
      m = message
      if condition(tocs):
         return tocs
      if not isinstance(m, str):
         m = m(tocs)
      raise pp.ParseException(s, loc, m)
  self.addParseAction(check_condition)
  return self

def addParseActionEx(self, pa, message = None):
  """
  Add parse action to a given pyparsing ParseElemenet,
  that, if it raise an exception, fail with a given message
  """

  def parse_action(s, loc, x):
      try:
        return pa(x)
      except Exception as e:
        if message:
           msg = '{message}\n{e}'
        else:
           msg = str(e)
        raise pp.ParseException(s, loc, msg).with_traceback(e.__traceback__) from e

  self.addParseAction(parse_action)

class White(pp.White):
  """ Fix for whitechars in pp.White

      In Python 3.10, pyparsing.White do not respect default_white_chars.
      This class fixes this.
  """

  def __init__(self, white):
      super().__init__(white)
      self.setWhitespaceChars('')

pp.ParserElement.addConditionEx = addConditionEx
pp.ParserElement.addParseActionEx = addParseActionEx
