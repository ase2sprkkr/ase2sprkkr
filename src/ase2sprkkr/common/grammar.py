"""
Various pyparsing grammar elements and a few useful routines.
"""

from contextlib import contextmanager
import pyparsing as pp
import re
from .decorators import cache
from typing import Optional


@contextmanager
def generate_grammar():
    """ Set the pyparsing newline handling
        and then restore to the original state """
    try:
      old = None
      kchars = None
      pe = pp.ParserElement
      if hasattr(pe, "DEFAULT_WHITE_CHARS"):
          old = pe.DEFAULT_WHITE_CHARS
      if hasattr(pp.Keyword, "DEFAULT_KEYWORD_CHARS"):
          kchars = pp.Keyword.DEFAULT_KEYWORD_CHARS + '-'
      pe.setDefaultWhitespaceChars(' \t\r')
      yield
    finally:
      if old is not None:
        pe.setDefaultWhitespaceChars(old)
      if kchars is not None:
        pp.Keyword.DEFAULT_KEYWORD_CHARS = kchars


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

  optional_quote = pp.Optional("'").suppress()
  """ Grammar for an optional quote """


def separator_pattern(char):
  return f'[{char}]' * 11 + '*'


@cache
def separator_grammar(char):
  """ Pattern for separating sections in an input file """
  separator = pp.Regex(separator_pattern(char)).setName(f"{char*10}[{char*4}....]").suppress()
  return separator


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
           msg = f'{message}\n{e}'
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


class SkipToRegex(pp.Token):
   """ Skip to given regex """

   def __init__(self, pattern, include_pattern:Optional[bool]=None,
                               parse_pattern:bool=False):
       if not isinstance(pattern, re.Pattern):
           pattern = re.compile(pattern)
       self.pattern = pattern
       if include_pattern is None:
           include_pattern = parse_pattern
       self.include_pattern = include_pattern
       self.parse_pattern = parse_pattern
       super().__init__()

   @property
   def custom_name(self):
       self.customName="skipTo{self.pattern}"

   def parseImpl(self, instr, loc, doActions = True):
       result = self.pattern.search(instr,loc)
       if result:
           start = result.start()
           out = instr[loc:start]
           if self.parse_pattern:
               out = (out, instr[start:result.end()])
           if self.include_pattern:
               loc = result.end()
           else:
               loc = result.start()
           return loc, pp.ParseResults(out)
       raise pp.ParseException(instr, loc, "Pattern {self.pattern} not found", self)


pp.ParserElement.addConditionEx = addConditionEx
pp.ParserElement.addParseActionEx = addParseActionEx
