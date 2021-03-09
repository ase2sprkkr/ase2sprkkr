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

class BaseGrammar:
   def grammar(self):
       """ Generate grammar with the correct settings of pyparsing """
       with generate_grammar():
         return self._grammar()

   def _tuple_with_my_name(self, expr, delimiter=None):
        """ Create the grammar returning tuple (self.name, <result of the expr>) """
        if self.name_in_grammar:
            name = pp.CaselessKeyword(self.name).setParseAction(lambda x: self.name)
            if delimiter:
              name += delimiter
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)
        out = name + expr
        out.setParseAction(lambda x: tuple(x))
        return out


with generate_grammar():
  optional_line_end = pp.Suppress(pp.LineEnd() | pp.WordStart() ).setName(' ')
  line_end = pp.Suppress(pp.LineEnd()).setName('\n')
  end_of_file = pp.Suppress(pp.ZeroOrMore(pp.LineEnd()) + pp.StringEnd()).setName('')
  separator = pp.Suppress(pp.Literal('*'*10) + pp.ZeroOrMore('*'))

def delimitedList(expr, delim):
  """ Delimited list with already suppressed delimiter (or with a in-results-wanted one) """
  return delim + pp.ZeroOrMore(delim + expr)
