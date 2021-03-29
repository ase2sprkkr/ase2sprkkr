import pyparsing as pp
import io
from pyparsing import Word, Suppress
import itertools
import functools
import numpy as np
from collections import Hashable
from .misc import OrderedDict
ppc = pp.pyparsing_common
from .grammar import generate_grammar, separator as separator_grammar, delimitedList, line_end
from .misc import class_property, copy_list

from ase.units import Rydberg
import copy

context =  generate_grammar()
context.__enter__()

class BaseType:
  """ Base class for definition of configuration option types """

  """ Common types have values, only the Separator does not """
  has_value = True

  """ Default value for BaseValueDefinition.name_in_grammar.
      Some types (e.g. Tables) commonly have no name (are identified
      by its position in the potential file.
  """

  name_in_grammar = True

  """ Default value for the given type. It can be overriden in the constructor. """
  default_value = None


  def __init__(self, prefix=None, postfix=None, format='', default_value=None, condition=None):
      self.prefix = prefix
      self.postfix = postfix
      self.format = format
      self.condition = condition

      """ Some subclasses has default_value defined via read-only property. """
      if default_value is not None:
         self.default_value = self.convert(default_value)

  def __str__(self):
      return self.__class__.__name__

  @functools.cached_property
  def grammar(self):
    """ Return a pyparsing grammar for the type """
    grammar = self._grammar
    if self.prefix:
       grammar = pp.Literal(self.prefix).suppress().setName('{') + grammar
    if self.postfix:
       grammar += pp.Literal(self.postfix).suppress().setName('}')
    grammar = self.transform_grammar(grammar)
    grammar.addCondition(lambda x: self.validate(x[0]))
    return grammar

  def grammar_name(self):
    """ Human readable expression of the grammar. By default,
        this is what is set by grammar.setName, however, sometimes
        is desirable to set even shorter string """
    return str(self.grammar)

  def transform_grammar(self, grammar):
    """ Chance for the resulting class to alter the resulting prefixed grammar """
    return grammar

  def missing_value(self):
    """ The =value flag can be ommited in the confing file
    Return
    ------
    can_be_ommited : bool
        Is an ommision of the value possible, e.g. the option is given as Flag (only by name of the option)
    default_value
        The value used if the value is ommitted
    do_not_output_value
        The value, with which the variable should not be outputed at all (e.g. False for a flag)
    """
    return False, None, None

  def validate(self, value, param_name='<Unknown>'):
    """ Validate either the pyparsing result or a user given value

    Paramters
    ---------
    value : mixed
      Value to be validated
    param_name : str
      Parameter name to be used in possible throwed exception (Optional)
    """
    err = self._validate(value)
    if err is not True:
      self._valueError(value, err, param_name)
    return True

  def _validate(self, value):
    """ Return error message if the value is not valid """
    return True

  def _valueError(self, value, error_message=False, param_name = '<Unkwnown>'):
    if error_message is False:
       error_message = 'invalid value'
    raise ValueError("Value '{}' for parameter {} is not valid: {}".format(value, param_name, error_message))

  def read(self, token, parameterName='<Unknown>'):
    """ Transform pyparsing token to a validated value """
    self.validate(val)
    return val

  def convert(self, value):
    """ Convert a value from user to a "cannonical form" """
    return value

  def _string(self, val):
    return str(val)

  def string(self, val):
    val = self._string(val)
    if self.prefix:
       val = self.prefix + val
    if self.postfix:
       val += self.postfix
    if self.print_width:
       spaces = self.print_width - len(val)
       if spaces > 0:
          val = " "*spaces + val
    return val

  def write(self, f, val):
    f.write(self.string(val))

  def print(self, val):
    s = io.StringIO()
    self.write(s, val)
    return s.getvalue()

  def copy(self):
    return copy.copy(self)

class Unsigned(BaseType):

  _grammar = ppc.integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(value, param_name):
    if not isinstance(value, int): return "Integer value required"
    return value >= 0 or "Positive value required"

  def grammar_name(self):
    return '<+int>'

  numpy_type = int

Unsigned.I = Unsigned()

class Integer(BaseType):

  _grammar = ppc.signed_integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(self, value):
    return isinstance(value, (int, np.int64) ) or f'Integer value required ({value.__class__} was given)'

  def grammar_name(self):
    return '<int>'

  numpy_type = int

Integer.I = Integer()

class Bool(BaseType):

  _grammar = (pp.Keyword('T') | pp.Keyword('F')).setParseAction( lambda x: x[0] == 'T' )

  def _validate(self, value):
    return isinstance(value, bool) or "Bool value rquired"

  def grammar_name(self):
    return '<T|F>'

  def _string(self, val):
    return 'T' if val else 'F'

  numpy_type = bool

Bool.I = Bool()


class Real(BaseType):

  _grammar = ppc.fnumber.setParseAction(lambda x: float(x[0]))

  def _validate(self, value):
    return isinstance(value, float) or "Float value required"

  def grammar_name(self):
    return '<float>'

  numpy_type = float

Real.I = Real()

class Date(BaseType):

  _grammar = pp.Regex(r'(\d{2}).(\d{2}).(\d{4})').setParseAction(lambda x: datetime.date(x[2], x[1], x[0]))

  def _validate(self, value):
    return isinstance(value, datetimei.date) or "Date (datetime.date) value required"

  def grammar_name(self):
    return '<dd.mm.yyyy>'

Date.I = Date()


class RealWithUnits(BaseType):

  @class_property
  def _grammar(cls):
    units = pp.Or((pp.Keyword(v).setParseAction(lambda x,*args, u=u: u)  for v,u in  cls.units.items()))
    units = units | pp.Empty().setParseAction(lambda x: 1.)
    cls._grammar = ppc.fnumber.setParseAction(lambda x: float(x[0])) + pp.Or(units)
    cls._grammar.setParseAction(lambda x: x[0]*x[1])
    return cls._grammar

  def _validate(self, value):
    return isinstance(value, float) or "Float value required"

  def grammar_name(self):
    return '<float>[{}]'.format("|".join(self.units))

  numpy_type = float

class Energy(RealWithUnits):

  units = {
      'Ry' : 1.,
      'eV' : 1. / Rydberg
  }

Energy.I = Energy()

class String(BaseType):

  _grammar = Word(pp.printables,excludeChars=",;{}").setParseAction(lambda x:x[0])

  def _validate(self, value):
    return isinstance(value, str) or "String value required"

  def grammar_name(self):
    return '<str>'

  numpy_type = object

String.I = String()

class QString(String):
  """ Quoted string"""
  _gramar = (pp.Word(pp.printables) or pp.QuotedString("'")).setParseAction(lambda x:x[0])

QString.I = QString()

class Keyword(BaseType):

  def __init__(self, *keywords, **kwargs):
    super().__init__(**kwargs)
    self.keywords = keywords
    self._grammar = pp.MatchFirst((pp.CaselessKeyword(i) for i in self.keywords)).setParseAction(lambda x: x[0].upper())

  def _validate(self, value):
    return value in self.keywords or "Required one of [" + "|".join(self.keywords) + "]"

  def grammar_name(self):
      return '|'.join(('"'+i+'"' for i in self.keywords ))

  def __str__(self):
      return self.grammar_name()

  def convert(self, value):
      return value.upper()

def DefKeyword(default, *others, **kwargs):
  return Keyword(default, *others, default_value=default, **kwargs)


class Flag(BaseType):

  def grammar_name(self):
      return None

  def str(self):
      return "(Flag)"

  def missing_value(self):
      return (True, True, False)

  def _validate(self, value):
      return value is True or value is False or value is None or "This is Flag with no value, please set to True to be present or to False/None to not"

  _grammar = pp.Empty().setParseAction(lambda x: True)

Flag.I = Flag()

normalize_type_map = {
    np.int64 : int,
    np.float64: float,
    np.bool: bool
}

def normalize_type(type):
    return normalize_type_map.get(type, type)

type_from_type_map = OrderedDict([
    (float, Real.I),
    (int  , Integer.I),
    (bool,  Bool.I),
    (str  , String.I)]
)

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

      Paramteres
      ----------
      type: A python type or BaseType
        A type to be converted to a grammar element (BaseType class descendant)

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

  delimiter=pp.White(' \t').suppress()
  delimiter_str = ' '

  def __init__(self, type, length=None, max_length=None, min_length=None, **kwargs):
    super().__init__(**kwargs)
    self.min_length = min_length or length
    self.max_length = max_length or length
    self.type = type_from_type(type)
    grammar = self.type.grammar()
    grammar = delimitedList(grammar, self.delimiter)
    grammar.setParseAction(lambda x: np.atleast_1d(x.asList()))
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

  def _validate(self, value):
    for i,v in enumerate(value):
        out = self.type._validate(v)
        if out is not True:
           return "Value {} in the set is incorrect: {}".format(i, out)
    if self.min_length is not None and len(value) < self.min_length:
       return f"The set shoud be at least {self.min_length} items long, it has {len(value)} items"
    if self.max_length is not None and len(value) > self.min_length:
       return f"The set can not have more than {self.max_length} items, it has {len(value)} items"
    return True

  def convert(self, value):
    if not isinstance(value, np.ndarray):
       return np.atleast_1d( value )
    return value


class SetOf(Array):
  """ Set of values, e.g. {1,2,3} """

  delimiter = pp.Suppress(pp.Literal(',') | pp.Literal(';') | pp.White(' \t')).setName('[,; ]')
  delimiter_str = ','

  def __init__(self, type, **kwargs):
    kwargs.setdefault('prefix', '{')
    kwargs.setdefault('postfix', '}')
    super().__init__(type, **kwargs)

  def transform_grammar(self, grammar):
    return grammar | self.type.grammar.copy().addParseAction(lambda x: np.atleast_1d(x.asList()))

  def __str__(self):
    return "SetOf({})".format(str(self.type))



""" Map python native types to configuration value types """

type_from_set_map = OrderedDict([
    (float, SetOf(float)),
    (int  , SetOf(int)),
])

def type_from_value(value):
  """ Gues the option type from a python value. E.g. 2 => Integer """
  if isinstance(value, (list, np.ndarray)):
     return type_from_set_map[normalize_type(value[0].__class__)] if len(value) else Integer.I
  if isinstance(value, str):
     try:
        String._grammar.parseString(value, True)
        return String.I
     except Exception:
        return QString.I

  return type_from_type(value.__class__).__class__(default_value = value)

def type_from_default_value(value, format='', format_all=False):
   if inspect.isclass(value) or isinstance(value, BaseType):
      return type_from_type(value, format=format, format_all=format_all)
   ptype = normalize_type(value.__class__)
   gtype = type_from_type(value.__class__).__class__
   return gtype(default_value = value, format=format_for_type(format, ptype))

class Mixed(BaseType):

  _grammar = pp.Or((
    i.grammar() for i in [
      Real.I,
      Integer.I,
      Bool.I,
      Energy.I,
      type_from_set_map[int],
      type_from_set_map[float],
      QString.I,
      String.I,
    ]
  ))

  def _validate(self, value):
    type = type_from_value(value)
    if type is value:
       return 'Can not determine the type of value {}'.format(value)
    return type.validate(value)

  def grammar_name(self):
    return '<mixed>'

  def _string(self, val):
    if isinstance(val, bool):
       return Bool._string(self, val)
    else:
       return super()._string(val)

Mixed.I = Mixed()

class Separator(BaseType):
  """ Special class for **** separator inside a section """

  _grammar = separator_grammar
  has_value = False

  def _validate(self, value):
      return 'Can not set a value to a separator'

  def _grammar_name(self):
      return '****...****\n'

  def _string(self, val=None):
      return '*'*80

Separator.I = Separator()

class Sequence(BaseType):
  """ A sequence of values of given types """

  def __init__(self, *types, format='', format_all=False, allowed_values=None, **kwargs):
      super().__init__(**kwargs)
      if isinstance(format, (str, dict)):
        format = itertools.repeat(format)
      self.types = [ type_from_default_value(i, dfs, format_all=format_all) for i,dfs in zip(types, format) ]
      self._grammar = pp.And([i.grammar for i in self.types]).setParseAction(lambda x: tuple(x))
      if allowed_values is not None:
         if not isinstance(allowed_values, set):
            allowed_values = set(allowed_values)
         self._grammar.addConditionEx(lambda x: x[0] in allowed_values, lambda x: f'{x[0]} is not in the list of allowed values')

  def _validate(self, value):
      if not isinstance(value, tuple) or len(value) != len(self.types):
          return f'A tuple of {len(self.types)} values is required'
      for i,j in zip(self.types, value):
          out = i.validate(j)
          if out is not True: return out
      return True

  def grammar_name(self):
      return  " ".join( (f'{j.grammar_name()}' for j in self.types) )

  def _string(self, val):
      out = []
      for i,v in zip(self.types, val):
          out.append(' ')
          out.append(i.string(v))
      return ''.join(out)


class Table(BaseType):
  """ Table, optionaly with named columns, e.g.
      IQ     IREFQ       IMQ       NOQ  ITOQ  CONC
       1         1         1         1     1 1.000
       2         2         2         1     2 1.000
  """

  name_in_grammar = False

  def __init__(self, columns=None, column_widths = 16, header=None,
                     numbering=None, numbering_label=None,
                     prefix=None, postfix=None, length=None,
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
      self.sequence = Sequence( *columns, column_widths = column_widths )
      self.header = header
      if numbering.__class__ is str:
         numbering_label=numbering
         numbering=True
      self.numbering = Integer.I if numbering is True else numbering
      self.numbering_label = numbering_label
      self.named_result = named_result

      line = self.sequence.grammar()
      if self.numbering:
         line = self.numbering.grammar() + line
      grammar = delimitedList(line, line_end)
      if self.names:
         grammar = pp.Suppress(pp.And(self.names) + pp.lineEnd) + grammar
         if self.numbering_label:
           grammar = pp.Keyword(self.numbering_label).suppress() + grammar

      def ensure_numbering(s, loc, x):
          numbers = x[::2]
          datas = x[1::2]
          if not numbers == [*range(1, len(numbers)+1)]:
             raise pp.ParseException(s, loc, 'First column should contain row numbering')
          return datas

      if self.numbering is not None:
         grammar.addParseAction(ensure_numbering)

      grammar.addParseActionEx( lambda x: np.array(x.asList(), self.numpy_type), "Cannot retype to numpy array")
      if length:
         grammar.addConditionEx(lambda x: len(x[0]) == length, lambda x: f'Just {length} rows are required, {len(x[0])} found.')
      self._grammar = grammar

  def _string(self, data):
      out = []
      if self.header:
         def gen():
             for n,t in zip(self.names, self.sequence.types):
                 yield n
                 yield t.print_width
         fstr = (" {:>{}}"*len(self.names))

         if self.numbering_label:
             fstr = String(print_width=self.numbering.print_width).string(self.numbering_label) + fstr
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


  @functools.cached_property
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

  def grammar_name(self):
      if self.names:
        data = " ".join( (f'{i}:{j.grammar_name()}' for i,j in zip(self.names, self.sequence.types) ) )
      else:
        data = self.sequence.grammar_name()
      return f"<TABLE of {data}>"

integer = Integer.I
unsigned = Unsigned.I
boolean = Bool.I
flag = Flag.I
real = Real.I
string = String.I
qstring = QString.I
mixed = Mixed.I
separator = Separator.I
energy = Energy.I


context.__exit__(None, None, None)
del context
