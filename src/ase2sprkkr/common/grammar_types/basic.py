""" Common GrammarTypes as numbers, strings etc. """

from ase.units import Rydberg
import datetime
import pyparsing as pp
from typing import Optional
import numpy as np
import re

from ..decorators import add_to_signature, cached_property
from ..grammar import generate_grammar, separator_grammar, \
                     replace_whitechars, optional_quote
from .grammar_type import TypedGrammarType, GrammarType, add_to_parent_validation

ppc = pp.pyparsing_common

class Char(TypedGrammarType):

  datatype = str

  @add_to_parent_validation
  def _validate(self, value, why='set'):
      return len(value) == 1 or "Char has to have length one"

  def grammar_name(self):
      return "<char>"

  _grammar = pp.Word(pp.printables, exact=1)


class Number(TypedGrammarType):
  """ Base class for a number - descendants of this class can have minimal and/or maximal possible value. """

  @add_to_signature(GrammarType.__init__)
  def __init__(self, min:Optional[int]=None, max:Optional[int]=None, *args, **kwargs):
      """
      Parameters
      ----------
      min:
        Minimal allowed value.

      max:
        Maximal allowed value.
      """
      self.min = min
      self.max = max
      super().__init__(*args, **kwargs)

  @add_to_parent_validation
  def _validate(self, value, why='set'):
      if self.min is not None and self.min > value:
         return f"A value greater that or equal to {self.min} is required, {value} have been given."
      if self.max is not None and self.max < value:
         return f"A value less than or equal to {self.max} is required, {value} have been given."
      return True


class FixedPointNumber(Number):

  numpy_type = int

  def convert(self, value):
      if isinstance(value, float) and np.abs(int(value) - value) < 1e-15:
          from .warnings import SuspiciousValueWarning
          SuspiciousValueWarning(value, 'An attemtp to set <integer> value using an <float>. I will do a conversion for you.').warn(
              stacklevel=6
          )
          value=int(value)

      return super().convert(value)

  @add_to_parent_validation
  def _validate(self, value, why='set'):
      if isinstance(value, float):
          return "A float value is not allowed for integer."
      return super()._validate(value, why)


class Unsigned(FixedPointNumber):
  """ Unsigned integer (zero is possible) """

  _grammar = replace_whitechars(ppc.integer).setParseAction(lambda x:int(x[0]))

  @add_to_parent_validation
  def _validate(self, value, why='set'):
      return value >= 0 or "A positive value required"

  def grammar_name(self):
    return '<+int>'

  datatype_name = 'unsigned integer'


class ObjectNumber(Unsigned):
  """ An abstract class, that describe an unsigned integer, that reffers to an object.
  User can give the object either using the object, or by the number. Descendant classes
  should take care of transforming the object to the resulting integer (by setting
  the result property of the described :class:`Option<ase2sprkkr.common.options.Option>`)

  The type of te object should be given by the type class property.
  """

  def convert(self, value):
      if isinstance(value, self.type):
         return value
      return super().convert(value)

  def _validate(self, value, why='set'):
      return isinstance(value, self.type) or super()._validate(value, why=why)


class Integer(FixedPointNumber):
  """ Signed integer """

  _grammar = replace_whitechars(ppc.signed_integer).setParseAction(lambda x:int(x[0]))

  def grammar_name(self):
    return '<int>'


class BaseBool(TypedGrammarType):
  """ A base type for all kind of boolean values."""
  numpy_type = bool
  type_name = 'boolean'


class Bool(BaseBool):
  """ A bool type, whose value is represented by a letter (T or F) """
  _grammar = (pp.CaselessKeyword('T') | pp.CaselessKeyword('F')).setParseAction( lambda x: x[0].upper() == 'T' )

  def grammar_name(self):
    return '<T|F>'

  def _string(self, val):
    return 'T' if val else 'F'


class Boolean(BaseBool):
  """ A bool type, whose value is represented by a letter (T or F) """
  _items = [ 'True', 'False', '1', '0', 'yes', 'no' ]
  _grammar =  pp.Or([pp.CaselessKeyword(i) for i in _items]).\
                   setParseAction( lambda x: x[0].lower() in ('true','yes','1') )

  def grammar_name(self):
    return '<True|False|0|1|yes|no>'

  def _string(self, val):
    return 'True' if val else 'False'


class IntBool(BaseBool):
  """ A bool type, whose value is represented by a letter (1 or 0) """
  _grammar = (pp.CaselessKeyword('1') | pp.CaselessKeyword('0')).setParseAction( lambda x: x[0] == '1' )
  _rev_grammar = _grammar.copy().setParseAction( lambda x: x[0] == '0' )

  @add_to_signature(TypedGrammarType.__init__)
  def __init__(self, reversed=False, *args, **kwargs):
      """
      Parameters
      ----------
      reversed
       "reversed integer-boolean" returns 1 if it is False
      """
      self.reversed = bool(reversed)
      super().__init__(*args, **kwargs)

  def grammar_name(self):
      return '<1|0>'

  def _string(self, val):
      return '1' if val != self.reversed else '0'


class Real(Number):
  """ A real value """
  _grammar = replace_whitechars(ppc.fnumber).setParseAction(lambda x: float(x[0]))

  def grammar_name(self):
    return '<float>'

  numpy_type = float

  nan = None

  @add_to_signature(Number.__init__)
  def __init__(self, nan=None ,*args, **kwargs):
      self.nan = nan
      if nan is not None:
          self._grammar = self._grammar | pp.Regex(nan).set_parse_action(lambda x: float('NaN'))
      super().__init__(*args, **kwargs)

  def convert(self, value):
      if isinstance(value, (int, np.integer)):
          from .warnings import SuspiciousValueWarning
          SuspiciousValueWarning(value, 'An attemtp to set <float> value using an <integer>. I will do a conversion for you.').warn(
              stacklevel=6
          )
          value=float(value)

      return super().convert(value)


class Date(Number):
  """ A date value of the form 'DD.MM.YYYY' """

  _grammar = pp.Regex(r'(?P<d>\d{2}).(?P<m>\d{2}).(?P<y>\d{4})').setParseAction(lambda x: datetime.date(int(x['y']), int(x['m']), int(x['d'])))

  def grammar_name(self):
    return '<dd.mm.yyyy>'

  def _string(self, val):
    return val.strftime("%d.%m.%Y")

  numpy_type = datetime.date
  type_name = 'date'


class BaseRealWithUnits(Real):
  """ The base class for float value, which can have units append.
      The value is converted automatically to the base units.
  """

  grammar_cache = {}
  """ The grammar for units is cached """

  def _grammar_units(self, units):
    i = id(units)
    if not i in self.grammar_cache:
      units = pp.Or(
        (pp.Empty() if v is None else pp.CaselessKeyword(v))
        .setParseAction(lambda x,*args, u=u: u) for v,u in units.items()
      )
      out = Real.I.grammar() + pp.Or(units)
      out.setParseAction(lambda x: x[0] * x[1])
      self.grammar_cache[i] = out
      return out
    return self.grammar_cache[i]

  def _grammar(self, param_name):
    return self._grammar_units(self.units)

  def _validate(self, value, why='set'):
    return isinstance(value, float) or "A float value required"

  def grammar_name(self):
    return '<float>[{}]'.format("|".join(('' if i is None else i for i in self.units)))

  numpy_type = float


class RealWithUnits(BaseRealWithUnits):
  """ A float value with user-defined units """

  def __init__(self, *args, units, **kwargs):
     self.units = units
     super().__init__(*args, **kwargs)


class Energy(BaseRealWithUnits):
  """ The grammar type for energy. The default units are Rydberg, one can specify eV. """

  units = {
      'Ry' : 1.,
      'eV' : 1. / Rydberg,
      None : 1.,
  }
  """ The allowed units and their conversion factors """

  def __str__(self):
      return "Energy (<Real> [Ry|eV])"


class BaseString(TypedGrammarType):
  """ Base type for string grammar types """

  datatype = str
  datatype_name = 'string'

  @add_to_parent_validation
  def _validate(self, value, why='set'):
    if not why=='parse':
      try:
        self._grammar.parseString(value, True)
      except pp.ParseException as e:
        return f"Forbidden character '{e.line[e.col-1]}' in the string"
    return True


class String(BaseString):
  """ Just a string (without whitespaces and few special chars) """
  _grammar = pp.Word(pp.printables,excludeChars=",;{}").setParseAction(lambda x:x[0])

  def grammar_name(self):
    return '<str>'

class AlwaysQString(BaseString):
  """ Either a quoted string, or just a word (without whitespaces or special chars). Always printed with quotes. """
  _grammar = (pp.Word(pp.printables, excludeChars="\",;{}") | pp.QuotedString('"')).setParseAction(lambda x:x[0])

  def _validate(self, value, why='set'):
    if not value.__class__ is str:
        return "String expected"
    return True

  def _string(self, value):
    value = value.replace('"', '\\"')
    return f'"{value}"'

  def grammar_name(self):
    return '"<str>"'


class QString(AlwaysQString):
  """ Either a quoted string, or just a word (without whitespaces or special chars) """

  def _string(self, value):
    value = str(value)
    if re.match(r"^[\S\"]+$", value):
        return value
    return super()._string(value)


class LineString(BaseString):
  """ A string, that takes all up to the end of the line """
  _grammar = pp.SkipTo(pp.LineEnd() | pp.StringEnd())

  def grammar_name(self):
    return "<str....>\n"


class Keyword(GrammarType):
  """
  A value, that can take values from the predefined set of strings.
  """
  def __init__(self, *keywords, aliases=None, transform='upper', quote=None, **kwargs):
    self.aliases = aliases or {}
    self.quote = quote
    if transform == 'upper':
        self.transform = lambda x: str(x).upper()
    elif transform == 'lower':
        self.transform = lambda x: str(x).lower()
    else:
        self.transform = lambda x:x

    if len(keywords)==1 and isinstance(keywords[0], dict):
       self.choices = keywords[0]
       keywords = self.choices.keys()
    else:
       self.choices = None

    self.keywords = [ self.transform(i) for i in keywords ]
    with generate_grammar():
      self._grammar = optional_quote + pp.MatchFirst((pp.CaselessKeyword(str(i)) for i in self.keywords)).setParseAction(lambda x: self.transform(x[0])) + optional_quote

    super().__init__(**kwargs)

  def _validate(self, value, why='set'):
    return value in self.keywords or "Required one of [" + "|".join(self.keywords) + "]"

  def _string(self, val):
      if self.quote and val.__class__ == str:
          return f"{self.quote}{val}{self.quote}"
      else:
          return str(val)

  def grammar_name(self):
      if len(self.keywords) == 1:
         return f'FixedValue({next(iter(self.keywords))})'
      return 'AnyOf(' + ','.join((str(i) for i in self.keywords )) + ')'

  def __str__(self):
      return self.grammar_name()

  def convert(self, value):
      out = self.transform(value)
      return self.aliases.get(out, out)

  def additional_description(self, prefix=''):
      ad = super().additional_description(prefix)
      if not self.choices:
         return ad
      out = f'\n{prefix}Possible values:\n'
      out += '\n'.join([f"{prefix}  {str(k):<10}{v}" for k,v in self.choices.items()])
      if ad:
         out += f'\n\n{prefix}' + ad
      return out

  is_independent_on_the_predecessor = True


def DefKeyword(default, *others, **kwargs):
  """
  A value, that can take values from the predefined set of strings, the first one is the default value.
  """
  if isinstance(default, dict) and len(others) == 0:
     def_val = next(iter(default))
  else:
     def_val = default
  return Keyword(default, *others, default_value=def_val, **kwargs)


class Flag(BaseBool):
  """
  A boolean value, which is True, if a name of the value appears in the input file.
  The resulting type of a Flag value is bool """

  _grammar = pp.Empty().setParseAction(lambda x: True)

  def grammar_name(self):
      return None

  def str(self):
      return "(Flag)"

  def missing_value(self):
      return (True, True, False)

  def _validate(self, value, why='set'):
      return value is True or value is False or value is None or "This is Flag with no value, please set to True to be present or to False/None to not"


class BasicSeparator(GrammarType):
  """ Basic type for separators - fake items in
  input/output file, which has no value """

  has_value = False

  def _validate(self, value, why='set'):
      return 'You can not set a value to a separator'

class KeywordSeparator(BasicSeparator):
  """ Separator with an exact string """
  def __init__(self, keyword, **kwargs):
      self.keyword = keyword
      super().__init__(**kwargs)

  @cached_property
  def _grammar(self):
      return pp.Keyword(self.keyword.strip())

  def _grammar_name(self):
      return self.keyword

  def _string(self, val=None):
      return self.keyword

class Separator(BasicSeparator):
  """ Special class for a separator inside a section.
      By default, it is a line of stars.
  """
  @cached_property
  def _grammar(self):
      return separator_grammar(self.char)

  def __init__(self, grammar=None, char='*', length=79, *args, **kwargs):
      self.char = char
      self.length = length
      if grammar:
         self._grammar = grammar
      super().__init__(*args, **kwargs)

  def _grammar_name(self):
      return f'{self.char*4}...{self.char*4}\n'

  def _string(self, val=None):
      return self.char * self.length


class BlankSeparator(BasicSeparator):
  """ Special class for a blank separator. In fact (with a delimiter) it is a blank line.
  """
  _grammar = pp.Empty().setParseAction(lambda x: [None])

  def _grammar_name(self):
      return ''

  def _string(self, val=None):
      return ''


# commonly used types
integer = Integer.I = Integer()
""" A standard grammar type instance for (signed) integers """
unsigned = Unsigned.I = Unsigned()         # NOQA: E741
""" A standard grammar type instance for unsigned integers """
boolean = Bool.I = Bool()                  # NOQA: E741
""" A standard grammar type instance for booleans in potential files """
flag = Flag.I = Flag()                     # NOQA: E741
""" A standard grammar type instance for booleans in input files """
real = Real.I = Real()                     # NOQA: E741
""" A standard grammar type instance for reals"""
date = Date.I = Date()                     # NOQA: E741
""" A standard instance for the grammar type for dates """
string = String.I = String()               # NOQA: E741
""" A standard grammar type instance for strings """
qstring = QString.I = QString()            # NOQA: E741
""" A standard grammar type instance for quoted strings in input files """
qstring = AlwaysQString.I = AlwaysQString() # NOQA: E741
""" A standard grammar type instance for quoted strings in configuration files """
line_string = LineString.I = LineString()  # NOQA: E741
""" A standard grammar type instance for one-line strings in potential files """
energy = Energy.I = Energy()               # NOQA: E741
""" A standard grammar type instance for energy values (float) for potential files """
separator = Separator.I = Separator()      # NOQA: E741
""" A standard grammar type instance for separators in potential files """
int_bool = IntBool.I = IntBool()           # NOQA: E741
""" A standard grammar type instance for bool expressed as integer """
boolean = Boolean.I = Boolean()            # NOQA: E741
""" A standard grammar type instance for True|False 0|1 yes|no boolean"""
