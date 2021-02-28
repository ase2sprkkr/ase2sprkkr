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
      pe.setDefaultWhitespaceChars(' \t')
      yield
    finally:
      if old is not None:
        pe.setDefaultWhitespaceChars(old)

class BaseGrammar:
   def grammar(self):
       with generate_grammar():
         return self._grammar()
