"""
Classes, that represents various value types that can appear in the configuration and problem definitionfiles.

Each grammar type can both parse string containing a value of a given type, and to create the string containing a given value.
"""


import pyparsing as pp
import io
import inspect
from pyparsing import Word, Suppress
import itertools
import numpy as np
from collections import namedtuple
from collections.abc import Hashable
from .misc import OrderedDict, cached_property, cache, add_to_signature
ppc = pp.pyparsing_common
from .grammar import generate_grammar, separator as separator_grammar, \
                     delimitedList, line_end, optional_quote,\
                     replace_whitechars, White

from ase.units import Rydberg
import copy
import datetime
from typing import Union, Any, Callable

context =  generate_grammar()
context.__enter__()
#it ensures that the generated grammar will have the correct whitespaces

class BaseType:
  """ Base class for definition of configuration option types

      A type without value (e.g. Separator) are just syntactical
      elements in the potentials file, that do not carry an information.
      Such elements do not yields (name, value) pair during parsing the file.

      Do not confuse this with BaseType.missing_value functionality.
      Missing_value is just the opposite: missing_value can be ommited in the file
      (or even the absence of the name in the file carry the information, that
      the Flag is False), but the name-value tuple of such Type is present
      in the parse result. On the other hand, has_value = False is in the file, but
      not in the result.

        """
  has_value = True

  """ Default value for BaseValueDefinition.name_in_grammar.
      Some types (e.g. Tables) commonly have no name (are identified
      by its position in the potential file.
  """
  name_in_grammar = True

  """ Default value for the given type. It can be overriden in the constructor (or just by setting
  the instantiated object attribute) """
  default_value = None

  """ Deafault type for creating numpy arrays (e.g. by Table) is object - to be redefined
      in the descendatns
  """
  numpy_type = object

  def __init__(self, prefix:Union[str,None]=None, postfix:Union[str,None]=None,
                     format:str='', default_value:Any=None,
                     condition:Union[Callable[[Any], Union[bool,str]],None]=None,
                     after_convert:Union[Callable[[Any], Any],None]=None):
      """
      Create the object.

      Parameters
      ----------
      prefix
        The string, that will be printed before the value

      postfix
        The string, that will be printed after the value

      format
        The (python) format string, that will be used for printing the value.
        The format is passed as format argument to ``str.format`` routine.

      default_value
        The default value of the options of this type. ``None`` means no default value.

      condition
        Function, that check the validity of the value. It should return ``True`` for a valid
        value, and ``False`` or string for invalid. The string is interpreted as an error message
        that explains the invalidity of the value.

      after_convert
        Function, that - if it is given - is applied to the (entered or parsed) value. The function
        is applied on the result of the
        :meth:`convert<ase2sprkkr.common.grammar_types.BaseType.convert>` method
      """

      self.prefix = prefix
      """ The string, that will be printed before the value """
      self.postfix = postfix
      """ The string, that will be printed after the value """
      self.format = format
      """ The (python) format string, that will be used for printing the value.
        The format is passed as format argument to ``str.format`` routine.  """
      self.condition = condition
      if after_convert is not None:
         self.convert = lambda v: \
              after_convert(self, self.__class__.convert(self, v))

      """ Some subclasses has default_value defined via read-only property. """
      if default_value is not None:
         self.default_value = self.convert(default_value)

  def __str__(self):
      return self.__class__.__name__

  @cache
  def grammar(self, param_name=False):
    """ Return a pyparsing grammar for the type """
    grammar = self._grammar
    if not isinstance(grammar, pp.ParserElement):
       grammar = grammar(param_name)
    if self.prefix or self.postfix:
       with generate_grammar():
        if self.prefix:
           grammar = pp.Literal(self.prefix).suppress().setName(self.prefix) + grammar
        if self.postfix:
           grammar += pp.Literal(self.postfix).suppress().setName(self.postfix)
        grammar = self.transform_grammar(grammar, param_name)

    if self.has_value:
       def validate(s, loc, x):
           try:
             out = self.validate(x[0], parse_check=True, param_name=param_name)
           except ValueError as e:
             raise pp.ParseException(s, loc, str(e) + '\nValidating of the parsed value failed') from e
           return x

       grammar.addParseAction(validate)
    grammar.grammar_type = self
    return grammar

  def parse(self, str, whole_string=True):
    return self.grammar().parseString(str, whole_string)[0]

  async def parse_from_stream(self, stream, up_to, start=None, whole_string=True):
    result = await stream.readuntil(up_to)
    result = result[:-len(up_to)].decode('utf8')
    if start:
       result = start + result
    return self.parse(result, whole_string)

  def grammar_name(self):
    """ Human readable expression of the grammar. By default,
        this is what is set by grammar.setName, however, sometimes
        is desirable to set even shorter string """
    return str(self.grammar)

  def transform_grammar(self, grammar, param_name=False):
    """ The chance for the resulting class to alter the resulting prefixed grammar """
    return grammar

  def missing_value(self):
    """ Is the configuraion value a flag? I.e. can be =<value> ommited
    in the configuration

    Return
    ------
    can_be_ommited : bool
        Is an ommision of the value possible, e.g. the option is given as Flag (only by name of the option)
    default_value
        The value used if the value is ommitted
    do_not_output_the_option
        The value, for which the variable should not be outputed at all (e.g. False for a flag)
    """
    return False, None, None

  def validate(self, value, param_name='<Unknown>', parse_check=False):
    """ Validate either the pyparsing result or a user given value

    Parameters
    ---------
    value : mixed
      Value to be validated
    param_name : str or callable
      Parameter name to be used in possible throwed exception (Optional)
      If it is callable, it should be a function that returns the param_name
    """
    try:
      err = self._validate(value, parse_check)
    except ValueError as err:
      self._valueError(value, err, param_name)
    if err is not True:
      self._valueError(value, err, param_name)
    if self.condition:
      err = self.condition(value)
      if err is not True:
        self._valueError(value, err, param_name)
    return True

  def _validate(self, value, parse_check=False):
    """ Return error message if the value is not valid """
    return True

  def _valueError(self, value, error_message=False, param_name=False):
    if callable(param_name):
       param_name = param_name()
    if param_name:
       param = f'for paramater {param_name} of type {self}'
    else:
       param = f'for type {self}'

    if error_message is False:
       error_message = 'invalid value'
    if isinstance(error_message, Exception):
      raise ValueError("Value '{}' {} is not valid: {}".format(value, param, error_message)) from error_message
    else:
      raise ValueError("Value '{}' {} is not valid: {}".format(value, param, error_message))

  def read(self, token, parameterName='<Unknown>'):
    """ Transform pyparsing token to a validated value """
    self.validate(val)
    return val

  def convert(self, value):
    """ Convert a value from user to a "cannonical form" """
    return value

  def _string(self, val):
    return val

  def string(self, val):
    val = self._string(val)
    if self.prefix:
       val = self.prefix + str(val)
    if self.postfix:
       val = str(val) + self.postfix
    if self.format:
       val = "{:{}}".format(val, self.format)
    return str(val)

  def write(self, f, val):
    f.write(self.string(val))

  def print(self, val):
    print(self.string(val))

  def copy(self):
    return copy.copy(self)

  def enrich(self, option):
    """ Some types can add properties to the options that have
    the type, e.g. see Sequence.enrich, which adds the ability to
    access the items of the sequence using [] """
    pass

  def __repr__(self):
    return "<{}>".format(self.__class__.__name__)

class Unsigned(BaseType):
  """ Unsigned integer (zero is possible) """

  _grammar = replace_whitechars(ppc.integer).setParseAction(lambda x:int(x[0]))

  def _validate(self, value, parse_check=False):
    if not isinstance(value, int): return "An integer value required"
    return value >= 0 or "A positive value required"

  def grammar_name(self):
    return '<+int>'

  numpy_type = int

Unsigned.I = Unsigned()

class Integer(BaseType):
  """ Signed integer """

  _grammar = replace_whitechars(ppc.signed_integer).setParseAction(lambda x:int(x[0]))

  def _validate(self, value, parse_check=False):
    return isinstance(value, (int, np.int64) ) or f'An integer value required, ({value.__class__} was given)'

  def grammar_name(self):
    return '<int>'

  numpy_type = int

Integer.I = Integer()

class Bool(BaseType):
  """ A bool type, whose value is represented by a letter (T or F) """
  _grammar = (pp.CaselessKeyword('T') | pp.CaselessKeyword('F')).setParseAction( lambda x: x[0] == 'T' )

  def _validate(self, value, parse_check=False):
    return isinstance(value, bool) or "A bool value rquired"

  def grammar_name(self):
    return '<T|F>'

  def _string(self, val):
    return 'T' if val else 'F'

  numpy_type = bool

Bool.I = Bool()


class Real(BaseType):
  """ A real value """
  _grammar = replace_whitechars(ppc.fnumber).setParseAction(lambda x: float(x[0]))

  def _validate(self, value, parse_check=False):
    return isinstance(value, float) or "A float value required"

  def grammar_name(self):
    return '<float>'

  numpy_type = float

Real.I = Real()

class Date(BaseType):
  """ A date value of the form 'DD.MM.YYYY' """

  _grammar = pp.Regex(r'(?P<d>\d{2}).(?P<m>\d{2}).(?P<y>\d{4})').setParseAction(lambda x: datetime.date(int(x['y']), int(x['m']), int(x['d'])))

  def _validate(self, value, parse_check=False):
    return isinstance(value, datetime.date) or "A date (datetime.date) value required"

  def grammar_name(self):
    return '<dd.mm.yyyy>'

  def _string(self, val):
    return val.strftime("%d.%m.%Y")

Date.I = Date()


class BaseRealWithUnits(BaseType):
  """ The base class for float value, which can have units append.
      The value is converted automatically to the base units.
  """

  grammar_cache = {}
  """ A grammar for units is cached """

  def _grammar_units(self, units):
    i = id(units)
    if not i in self.grammar_cache:
      units = pp.Or(
        (pp.Empty() if v is None else pp.CaselessKeyword(v))
                .setParseAction(lambda x,*args, u=u: u) for v,u in  units.items()
        )
      out =  Real.I.grammar() + pp.Or(units)
      out.setParseAction(lambda x: x[0]*x[1])
      self.grammar_cache[i] = out
      return out
    return self.grammar_cache[i]

  def _grammar(self, param_name):
    return self._grammar_units(self.units)

  def _validate(self, value, parse_check=False):
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

Energy.I = Energy()

class BaseString(BaseType):
  """ Base type for string grammar types """

  def _validate(self, value, parse_check=False):
    if not isinstance(value, str): return "A string value required"
    if not parse_check:
      try:
        self._grammar.parseString(value, True)
      except pp.ParseException as e:
        return f"Forbidden character '{e.line[e.col]}' in the string"
    return True

class String(BaseString):
  """ Just a string (without whitespaces and few special chars) """
  _grammar = Word(pp.printables,excludeChars=",;{}").setParseAction(lambda x:x[0])

  def grammar_name(self):
    return '<str>'
String.I = String()

class QString(BaseString):
  """ Either a quoted string, or just a word (without whitespaces or special chars) """
  _grammar = (pp.Word(pp.printables, excludeChars=",;{}") or pp.QuotedString("'")).setParseAction(lambda x:x[0])

  def grammar_name(self):
    return "'<str>'"

QString.I = String()

class LineString(BaseString):
  """ A string, that takes all up to the end of the line """
  _grammar = pp.SkipTo(pp.LineEnd() | pp.StringEnd())

  def grammar_name(self):
    return "'<str....>\n'"

LineString.I = LineString()

class Keyword(BaseType):
  """
  A value, that can take values from the predefined set of strings.
  """

  def __init__(self, *keywords, **kwargs):
    super().__init__(**kwargs)
    self.keywords = [ i.upper() for i in keywords ]
    with generate_grammar():
      self._grammar = optional_quote + pp.MatchFirst((pp.CaselessKeyword(i) for i in self.keywords)).setParseAction(lambda x: x[0].upper()) + optional_quote

  def _validate(self, value, parse_check=False):
    return value in self.keywords or "Required one of [" + "|".join(self.keywords) + "]"

  def grammar_name(self):
      return '|'.join(('"'+i+'"' for i in self.keywords ))

  def __str__(self):
      return self.grammar_name()

  def convert(self, value):
      return value.upper()

def DefKeyword(default, *others, **kwargs):
  """
  A value, that can take values from the predefined set of strings, the first one is the default value.
  """
  return Keyword(default, *others, default_value=default, **kwargs)


class Flag(BaseType):
  """
  A boolean value, which is True, if a name of the value appears in the input file.
  """

  def grammar_name(self):
      return None

  def str(self):
      return "(Flag)"

  def missing_value(self):
      return (True, True, False)

  def _validate(self, value, parse_check=False):
      return value is True or value is False or value is None or "This is Flag with no value, please set to True to be present or to False/None to not"

  _grammar = pp.Empty().setParseAction(lambda x: True)

Flag.I = Flag()

normalize_type_map = {
    np.int64 : int,
    np.float64: float,
    np.bool_: bool
}
""" Mapping of alternative types to the 'canonical ones'. """

def normalize_type(type):
    """ Return the 'canonical type' for a given type

    I.e. it maps numpy internal types to standard python ones

    doctest:
    >>> normalize_type(np.int64)
    <class 'int'>
    """
    return normalize_type_map.get(type, type)

type_from_type_map = OrderedDict([
    (float, Real.I),
    (int  , Integer.I),
    (bool,  Bool.I),
    (str  , String.I)]
)
""" The standard grammar_types for python types.

The value type can be given by a standard python type, this map maps the
python type for the appropriate grammar_type class.
"""

def format_for_type(format, type):
  """
  Returns the format appropriate to the given type

  Parameters
  ----------
  format: str or dict
    If it is str, just return it.
    Dict should has the form { type : format_for_the_type } + { None : default_format }
  """
  if isinstance(format, dict):
     if type in format:
        return format[type]
     return format[None]
  return format

def type_from_type(type, format='', format_all=False):
  """ Guess and return the grammar element (BaseType class descendatnt) from a python type. E.g. int => Integer.

      The given format can be optionally set to the returned grammar element.

      Parameters
      ----------
      type: A python type or BaseType
        A type to be converted to a grammar type (BaseType class descendant)

      format: str or dict
        The format to be applied to the resulting class. If dict is given, see 'format_for_type'
        for the way how the format is determined

      format_all: boolean
        If False (default), the format is not applied, if instance of BaseType is given as
        the type parameter. Otherwise, a copy of the input type with the applied format is returned
  """
  if isinstance(type, Hashable) and type in type_from_type_map:
    type = normalize_type(type)
    format = format_for_type(format, type)
    type = type_from_type_map[type]
    if format:
        type = type.copy()
        type.format = format
    return type
  elif format_all:
    type = type.copy()
    type.format = format_for_type(format, normalize_type(type.numpy_type))
  return type


class Array(BaseType):
  """ A (numpy) array of values of one type """

  delimiter=White(' \t').suppress()
  delimiter_str = ' '

  def __init__(self, type, default_value=None,
               length=None, max_length=None, min_length=None,
               as_list=False,
               **kwargs):
    """
    Parameters
    ----------
    type
      The grammar type of the values in the list (it can be given by a python type)

    default_value
      The default value for the list

    length
      If it is set, the list have to have just this length (it sets ``min_`` and ``max_length`` to the ``length``)

    min_length
      The minimal allowed length of the list.

    max_length
      The maximal allowed length of the list.

    as_list
      Type of the value array. True means List, False means np.ndarray, or custom type (e.g. tuple)
      can be provided. However, the value can be set using tuple or list anyway.
    """
    if isinstance(type, (list, np.ndarray)):
        if default_value is not None:
           raise ValueException("It is not possible for an Array to provide default_value both in 'default_value' and in 'type' argument")
        default_value = type
        type = type[0].__class__
    self.type = type_from_type(type)
    self.as_list = as_list
    super().__init__(default_value=default_value, **kwargs)
    self.min_length = min_length or length
    self.max_length = max_length or length
    with generate_grammar():
      grammar = self.type.grammar()
      grammar = delimitedList(grammar, self.delimiter)
      if self.as_list:
        if callable(self.as_list):
          grammar = grammar.setParseAction(lambda x: self.as_list(x.asList()))
        else:
          grammar = grammar.setParseAction(lambda x: [x.asList()])
      else:
        grammar.setParseAction(lambda x: self.convert(x.asList()))
        grammar.setName(self.grammar_name())
    self._grammar = grammar

  def __str__(self):
    return "Array({})".format(str(self.type))

  def grammar_name(self):
      gn = self.type.grammar_name()
      if self.min_length is not None and self.min_length == self.max_length:
        return f'{self.min_length}*{gn}'
      return f'{gn}{self.delimiter_str}{gn}{self.delimiter_str}...'

  def _string(self, val):
    it = iter(val)
    i = next(it)
    out = self.type.string(i)
    for i in it:
       out += self.delimiter_str
       out += self.type.string(i)
    return out

  def _validate(self, value, parse_check=False):
    if callable(self.as_list):
       cls = self.as_list
    elif self.as_list:
       cls = list
    else:
       cls = np.ndarray
    if not isinstance(value, cls):
       return f'A value of the {cls} type is required, a {value.__class__} is given'

    for i,v in enumerate(value):
        try:
          self.type.validate(v, parse_check=False)
        except ValueError as e:
           raise ValueError("Value {} in the set is incorrect: {}".format(i, str(e))) from e
    if self.min_length is not None and len(value) < self.min_length:
       return f"The array should be at least {self.min_length} items long, it has {len(value)} items"
    if self.max_length is not None and len(value) > self.max_length:
       return f"The array can not have more than {self.max_length} items, it has {len(value)} items"
    return True

  def convert(self, value):
    if self.as_list:
       if callable(self.as_list):
          return value if isinstance(value, self.as_list) else self.as_list(value)
       else:
          return list(value) if isinstance(value, tuple) else value
    if not isinstance(value, np.ndarray):
       if self.type.numpy_type == object:
          #https://stackoverflow.com/questions/60939396/forcing-a-creation-of-1d-numpy-array-from-a-list-array-of-possibly-iterable-obje
          out = np.empty(len(value), object)
          out[:] = value
          return out
       else:
          return np.atleast_1d(value)
    return value


class SetOf(Array):
  """ Set of values of the same type. E.g. {1,2,3} """

  delimiter = pp.Suppress(pp.Literal(',') | pp.Literal(';') | White(' \t')).setName('[,; ]')
  delimiter_str = ','

  def __init__(self, type, **kwargs):
    kwargs.setdefault('prefix', '{')
    kwargs.setdefault('postfix', '}')
    super().__init__(type, **kwargs)

  def transform_grammar(self, grammar, param_name=False):
    return grammar | self.type.grammar(param_name).copy().addParseAction(lambda x: np.atleast_1d(x.asList()))

  def __str__(self):
    return "SetOf({})".format(str(self.type))

type_from_set_map = OrderedDict([
    (float, SetOf(float)),
    (int  , SetOf(int)),
])
""" Map the python type of a collection member to a grammar type of the collection.

Only canonical types are expected, see :meth:`ase2sprkkr.common.grammar_types.normalize_type`
"""

recognized_set_types = ( list, tuple, np.ndarray )

def type_from_value(value):
  """ Gues the grammar type from a python value.

  ..doctest::
  >>> type_from_value(2)
  <Integer>

  >>> type_from_value(2.0)
  <Real>
  """
  if isinstance(value, recognized_set_types):
     return type_from_set_map[normalize_type(value[0].__class__)] if len(value) else Integer.I
  if isinstance(value, str):
     try:
        String._grammar.parseString(value, True)
        return String.I
     except Exception:
        return QString.I
  type = type_from_type(value.__class__)
  if type is value.__class__:
     raise ValueError('Cannot determine grammar type from value {value}')
  return type.__class__(default_value = value)

def type_from_default_value(value, format='', format_all=False):
   """ Guess the grammar type from a value, that will become the default value of the grammar type.

   It has to create a new object instance, as it has to set the default
   value property of the returned object. An (output) format can be applied to the
   resulting grammar type

   Grammar types passed as types are left as is, unless format_all flag is set.
   """
   if inspect.isclass(value) or isinstance(value, BaseType):
      return type_from_type(value, format=format, format_all=format_all)
   ptype = normalize_type(value.__class__)
   gtype = type_from_type(value.__class__).__class__
   return gtype(default_value = value, format=format_for_type(format, ptype))

class BaseMixed(BaseType):
  """
  A variant type - it can hold "anything".
  """

  type = None
  """ The types, that the value can hold. To be redefined in the descendants. """

  string_type = None
  """ Type of string grammar_type to be used.  To be redefined in the descendants. """

  def _grammar(self, param_name=False):
      return pp.MatchFirst((
        i.grammar(param_name) for i in self.types
      ))

  def get_type(self, value):
      """ Return the type of the value """
      return self.string_type if isinstance(value, str) else type_from_value(value)

  def _validate(self, value, parse_check=False):
      type = self.get_type(value)
      if type is value:
         return 'Can not determine the type of value {}'.format(value)
      return type.validate(value, parse_check)

  def grammar_name(self):
      return '<mixed>'

  def convert(self, value):
      return self.get_type(value).convert(value)

class Range(BaseMixed):
  """ A range type - it accepts either one value or range of two values of a given type."""

  @add_to_signature(BaseMixed.__init__, prepend=True)
  def __init__(self, type, *args, **kwargs):
      self._type = type_from_type(type)
      super().__init__(*args, **kwargs)

  @cached_property
  def types(self):
      return [
          self._type,
          SetOf(self._type, min_length=2, max_length=2)
      ]

  def get_type(self, value):
      return self.types[1 if isinstance(value, recognized_set_types) else 0]

class Mixed(BaseMixed):
  """ A variant value to be used in input files (in unknown - custom - options) """

  string_type = QString.I
  """ Input files use quoted strings. """

  types = [
      Energy.I,
      Real.I,
      Integer.I,
      type_from_set_map[int],
      type_from_set_map[float],
      QString.I,
      Flag.I,
  ]

  def missing_value(self):
    return True, True, False

Mixed.I = Mixed()

class PotMixed(BaseMixed):
  """ A variant value to be used in potential files (in unknown - custom - options) """

  string_type = LineString.I
  """ Potential files use line strings. """

  types = [
      Energy.I,
      Real.I,
      Integer.I,
      Bool.I,
      type_from_set_map[int],
      type_from_set_map[float],
      LineString.I,
  ]

  def _string(self, val):
    if isinstance(val, bool):
       return Bool._string(self, val)
    else:
       return super()._string(val)

PotMixed.I = PotMixed()


class Separator(BaseType):
  """ Special class for ``****`` separator inside a section """

  _grammar = separator_grammar.copy().setParseAction(lambda x: [None])
  has_value = False

  def _validate(self, value, parse_check=False):
      return 'Can not set a value to a separator'

  def _grammar_name(self):
      return '****...****\n'

  def _string(self, val=None):
      return '*'*79

Separator.I = Separator()

class Sequence(BaseType):
  """ A sequence of values of given types """

  def __init__(self, *types, format='', format_all=False, allowed_values=None,
               default_values=False, names=None, **kwargs):
      super().__init__(**kwargs)
      if names:
         self.names = names if isinstance(names, dict) else {name:i for i,name in enumerate(names)}
         self.value_type = namedtuple("_".join(names), names)
         self.value_constructor = self.value_type
      else:
         self.names = None
         self.value_type = tuple
         self.value_constructor = lambda *x: tuple(x)
      if isinstance(format, (str, dict)):
        format = itertools.repeat(format)
      self.types = [ type_from_default_value(i, dfs, format_all=format_all) for i,dfs in zip(types, format) ]
      if allowed_values and not isinstance(allowed_values, set):
         allowed_values = set(allowed_values)
      self.allowed_values=allowed_values
      self.default_values=default_values

  def _grammar(self, param_name = False):
      def grm(type):
          g = type.grammar(param_name)
          if self.default_values and type.default_value is not None:
             g = g | pp.Empty().setParseAction(lambda x: type.default_value)
          return g

      grammars = [grm(i) for i in self.types]
      grammar = pp.And(grammars).setParseAction(lambda x: self.value_constructor(*x))
      if self.allowed_values is not None:
         grammar.addConditionEx(lambda x: x[0] in self.allowed_values, lambda x: f'{x[0]} is not in the list of allowed values')
      return grammar

  def _validate(self, value, parse_check=False):
      if not isinstance(value, (self.value_type)) or len(value) != len(self.types):
          return f'A tuple of {len(self.types)} values is required'
      for i,j in zip(self.types, value):
          out = i.validate(j, parse_check=parse_check)
      return True

  def convert(self, value):
      if not isinstance(value, self.value_type):
         return self.value_constructor(*value)
         try:
            return self.value_constructor(*value)
         except TypeError:
            pass
      return value

  def grammar_name(self):
      return  " ".join( (f'{j.grammar_name()}' for j in self.types) )

  def _string(self, val):
      out = []
      for i,v in zip(self.types, val):
          out.append(' ')
          out.append(i.string(v))
      return ''.join(out)

  def enrich(self, option):

      class cls(option.__class__):
         def _get_index(sulf, name):
           if self.names and isinstance(name, str):
              return self.names[name]
           return name

         def __getitem__(self, key):
           key = self._get_index(key)
           return self()[key]

         def __setitem__(self, key, value):
           key = self._get_index(key)
           v = list(self())
           v[key] = value
           self.set(v)

      if self.names:
        for n,i in self.names.items():
            (lambda i: setattr(cls, n, property(
                lambda self: self[i],
                lambda self, v: self.__setitem__(i, v)
            )))(i)

      option.__class__ = cls


class Table(BaseType):
  """
  Table, optionaly with named columns, e.g.

    ::text

      IQ     IREFQ       IMQ       NOQ  ITOQ  CONC
       1         1         1         1     1 1.000
       2         2         2         1     2 1.000

  """

  name_in_grammar = False

  def __init__(self, columns=None,
                     header=None, free_header=False,
                     format = {float: '>21.17', None: '>16'}, format_all=True,
                     numbering=None, numbering_label=None, numbering_format=True,
                     prefix=None, postfix=None, length=None,
                     row_condition=None,
                     default_values=False,
                     named_result = False, **kwargs):
      super().__init__(prefix=None, postfix=None)
      if columns is None:
         columns = kwargs
      if isinstance(columns, dict):
         self.names = list(columns.keys())
         columns = columns.values()
      else:
         self.names = None
      if header is None:
         header = self.names
      self.sequence = Sequence( *columns, format=format, format_all=format_all, condition = row_condition, default_values=default_values )
      self.header = header
      self.free_header = free_header
      if numbering.__class__ is str:
         numbering_label=numbering
         numbering=True
      self.numbering = Unsigned.I if numbering is True else numbering
      if self.numbering and numbering_format and not (numbering_format is True and self.numbering.format):
         if numbering_format is True:
            numbering_format = '<4'
         self.numbering = self.numbering.copy()
         self.numbering.format = numbering_format
      self.numbering_label = numbering_label
      self.named_result = named_result
      self.length = length

  def _grammar(self, param_name=False):
      line = self.sequence.grammar(param_name)
      if self.numbering:
         line = self.numbering.grammar() + line # + pp.And._ErrorStop()
      grammar = delimitedList(line, line_end)
      if self.names:
         if self.free_header:
             fh = pp.SkipTo(line_end) + line_end
             if callable(self.free_header):
               fh.addConditionEx(lambda x: self.free_header(x[0]),
                                    lambda x: f"This is not an allowed header for table {param_name}: {x[0]}" )
             grammar = pp.Suppress(fh) + grammar
         else:
             def names():
                for n in self.names:
                    if ' ' in n:
                      """ multiple column headers for one column are allowed
                          -- see Occupation section"""
                      yield from map(pp.CaselessKeyword, n.split(' '))
                    else:
                      yield pp.CaselessKeyword(n)

             grammar = pp.Suppress(pp.And(list(names())) + pp.lineEnd) + grammar
             if self.numbering_label:
               grammar = pp.CaselessKeyword(self.numbering_label).suppress() + grammar

      def ensure_numbering(s, loc, x):
          numbers = x[::2]
          datas = x[1::2]
          if not numbers == [*range(1, len(numbers)+1)]:
             raise pp.ParseException(s, loc, 'First column should contain row numbering')
          return datas

      if self.numbering is not None:
         grammar.addParseAction(ensure_numbering)

      grammar.addParseActionEx( lambda x: np.array(x.asList(), self.numpy_type), "Cannot retype to numpy array")
      return grammar

  def _string(self, data):
      out = []
      if self.header:
         def gen():
             names = ((i[1] if isinstance(i, tuple) else i) for i in self.names)

             for n,t in zip(self.names, self.sequence.types):
                 yield n
                 yield t.format
         fstr = (" {:{}}"*len(self.names))

         if self.numbering:
            fstr = self.numbering.string(self.numbering_label or '') + fstr
         else:
             fstr = fstr[1:]
         out.append(fstr.format(*gen()))
         newline = True
      else:
         newline = False

      line = 1
      for i in data:
         if newline:
            out.append('\n')
         newline = True
         if self.numbering is not None:
            out.append(self.numbering.string(line))
            line+=1
         out.append(self.sequence.string(i))
      return ''.join(out)

  def _validate(self, value, parse_check=False):
      if not isinstance(value, np.ndarray):
         return f"Numpy array as a value required {value.__class__} given"
      dtype = self.numpy_type
      dim = 1 if isinstance(dtype, list) else 2
      if len(value.shape) != dim:
         return f"The array should have dimension={dim}, it has dimension {len(value.shape)}"
      if value.dtype != self.numpy_type:
         return f"The data type of the value should be {dtype}, it is {value.dtype}"
      if dim==2 and value.shape[1] != len(self.sequence.types):
         return f"The array is required to have {len(self.sequence.types)} columns, it has {value.shape[1]}"
      if self.length is not None and self.length != value.shape[0]:
         return f"The array is required to have {self.length} rows, it has {value.shape[1]}"
      return True

  def convert(self, value):
      return np.asarray(value, dtype = self.numpy_type)

  @cached_property
  def numpy_type(self):
      types = self.sequence.types
      nr = self.names and self.named_result
      if not nr:
         dtype = types[0].numpy_type
         for t in types[1:]:
             if t.numpy_type != dtype:
                 nr = True
                 break
         else:
             return dtype
      names = self.names or itertools.repeat('')
      return list(zip(names, (i.numpy_type for i in types)))

  def number_of_collumns(self):
      return len(self.sequence.types)

  def zero_data(self, length):
      """ Return array of zeros with the given number of rows and
          with the dtype of the table
      """
      dtype = self.numpy_type
      if isinstance(dtype, list):
         return np.zeros(length, dtype)
      else:
         return np.zeros((length, self.number_of_collumns()), dtype)

  def grammar_name(self):
      if self.names:
        data = " ".join( (f'{i}:{j.grammar_name()}' for i,j in zip(self.names, self.sequence.types) ) )
      else:
        data = self.sequence.grammar_name()
      return f"<TABLE of {data}>"

integer = Integer.I
""" A standard signed integer grammar type instance """
unsigned = Unsigned.I
""" A standard unsigned integer grammar type instance """
boolean = Bool.I
""" A standard bool grammar type instance (for potential files) """
flag = Flag.I
""" A standard bool grammar type instance (for input files) """
real = Real.I
""" A standard real grammar type instance """
string = String.I
""" A standard string grammar type instance """
qstring = QString.I
""" A standard quoted string grammar type instance (for input files) """
line_string = LineString.I
""" A standard line string grammar type instance (for potential files) """
mixed = Mixed.I
""" A standard variant grammar type instance (for input files) """
pot_mixed = PotMixed.I
""" A standard variant grammar type instance (for potential files) """
separator = Separator.I
""" A standard separator line grammar type instance (for potential files) """
energy = Energy.I
""" A standard energy float value type instance (for potential files) """

context.__exit__(None, None, None)
del context
