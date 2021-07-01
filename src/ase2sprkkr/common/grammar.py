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
  line_end = pp.Suppress(pp.LineEnd()).setName('\n')
  end_of_file = (pp.Regex('[\s]*') + pp.StringEnd()).suppress().setName('<EOF>')
  pp.Suppress(pp.ZeroOrMore(pp.LineEnd()) + pp.StringEnd()).setName('')
  separator_pattern = r'\*'*10+'\**'
  separator = pp.Regex(separator_pattern).setName("**********[***....]").suppress()
  separator.pattern = separator_pattern
  optional_quote = pp.Optional("'").suppress()

def delimitedList(expr, delim):
  """ Delimited list with already suppressed delimiter (or with a in-results-wanted one) """
  return expr + pp.ZeroOrMore(delim + expr)

def addConditionEx(self, condition, message):
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
pp.ParserElement.addConditionEx = addConditionEx
pp.ParserElement.addParseActionEx = addParseActionEx
