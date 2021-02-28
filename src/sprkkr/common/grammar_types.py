import pyparsing as pp
import io
from pyparsing import Word, Suppress, delimitedList
import itertools
import functools
import numpy as np
from collections import OrderedDict
ppc = pp.pyparsing_common

class Base:
  """ Base class for definition of configuration option types """

  def __str__(self):
    return self.__class__.__name__

  @functools.cached_property
  def grammar(self):
    """ Return a pyparsing grammar for the type """
    self._grammar.addCondition(lambda x: self.validate(x[0]))
    return self._grammar

  def is_optional(self):
    """ The =value flag can be ommited in the confing file"""
    return False

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

  def write(self, f, val):
    f.write(str(val))

  def print(self, val):
    s = io.StringIO()
    self.write(s, val)
    return s.getvalue()

class Unsigned(Base):

  _grammar = ppc.integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(value, param_name):
    if not isinstance(value, int): return "Integer value required"
    return value >= 0 or "Positive value required"

  def grammar_name(self):
    return '<+int>'

  numpy_type = int

Unsigned.I = Unsigned()

class Integer(Base):

  _grammar = ppc.signed_integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(self, value):
    return isinstance(value, (int, np.int64) ) or f'Integer value required ({value.__class__} was given)'

  def grammar_name(self):
    return '<int>'

  numpy_type = int

Integer.I = Integer()

class Bool(Base):

  _grammar = (pp.Keyword('T') | pp.Keyword('F')).setParseAction( lambda x: x[0] == 'T' )

  def _validate(self, value):
    return isinstance(value, bool) or "Bool value rquired"

  def grammar_name(self):
    return '<T|F>'

  numpy_type = bool

Bool.I = Bool()


class Real(Base):

  _grammar = ppc.fnumber.setParseAction(lambda x: float(x[0]))

  def _validate(self, value):
    return isinstance(value, float) or "Float value required"

  def grammar_name(self):
    return '<float>'

  numpy_type = float

Real.I = Real()

class String(Base):

  _grammar = Word(pp.printables).setParseAction(lambda x:x[0])

  def _validate(self, value):
    return isinstance(value, str) or "String value required"

  def grammar_name(self):
    return '<str>'

  numpy_type = str

String.I = String()

class QString(String):
  """ Quoted string"""
  _gramar = (pp.Word(pp.printables) or pp.QuotedString("'")).setParseAction(lambda x:x[0])

QString.I = QString()

class Keyword(Base):

  def __init__(self, *keywords):
    self.keywords = keywords
    self._grammar = pp.MatchFirst((pp.CaselessKeyword(i) for i in self.keywords)).setParseAction(lambda x: x[0].upper())

  def _validate(self, value):
    return value in keywords or "Required one of [" + "|".join(self.keywords) + "]"

  def grammar_name(self):
      return '|'.join(('"'+i+'"' for i in self.keywords ))

  def __str__(self):
      return self.grammar_name()

class Flag(Base):

  def grammar_name(self):
      return None

  def str(self):
      return "(Flag)"

  def is_optional(self):
      return True

  def _validate(self, value):
      return value is True or value is False or value is None or "This is Flag with no value, please set to True to be present or to False/None to not"

  _grammar = pp.Empty().setParseAction(lambda x: True)
      
Flag.I = Flag()


type_from_type_map = OrderedDict([
    (int  , Integer.I),
    (np.int64  , Integer.I),
    (float, Real.I),
    (np.float64, Real.I),
    (bool,  Bool.I),
    (np.bool_,  Bool.I),
    (str  , String.I)]
)

def type_from_type(type):
  """ Gues the option type from a python type. E.g. int => Integer """
  if type in type_from_type_map:
    return type_from_type_map[type]
  return type



class SetOf(Base):
  """ Set of values, e.g. {1,2,3} """
  def __init__(self, type, length=None, max_length=None):
    self.min_length = length
    self.max_length = max_length or length
    self.type = type_from_type(type)
    self._grammar = (Suppress("{") + delimitedList(self.type.grammar()) + Suppress("}")).setParseAction(lambda x: np.array(x.asList()))

  def __str__(self):
    return "SetOf({})".format(str(self.type))

  def grammar_name(self):
      return '{' + self.type.grammar_name() + ",...}"

  def write(self, f,val):
    f.write('{')
    if not isinstance(val, (list, tuple)):
      val = [ val ]
    it = iter(val)
    v = next(it)
    self.type.write(f, v)
    for i in it:
       f.write(",")
       self.type.write(f, val)
    f.write('}')

  def _validate(self, value):
    for i,v in enumerate(value):
        out = self.type._validate(v)
        if out is not True:
           return "Value {} in the set is incorrect: {}".format(i, out)
    return True

  def convert(self, value):
    if self.length is None and not isinstance(list(value)):
      return [ value ]


""" Map python native types to configuration value types """


type_from_set_map = OrderedDict([
    (int  , SetOf(int)),
    (float, SetOf(float))
])
type_from_set_map[np.int64] = type_from_set_map[int]
type_from_set_map[np.float64] = type_from_set_map[float]


def type_from_value(value):
  """ Gues the option type from a python value. E.g. 2 => Integer """
  if isinstance(value, (list, np.ndarray)):
     return type_from_set_map[value[0].__class__] if len(value) else Integer.I
  return type_from_type(value.__class__)

class Mixed(Base):

  _grammar = pp.Or((
    i.grammar() for i in
    itertools.chain(
      type_from_set_map.values(),
      type_from_type_map.values(), 
    )
  ))

  def _validate(self, value):
    type = type_from_value(value)
    if type is value:
       return 'Can not determine the type of value {}'.format(value)
    return type.validate(value)

  def grammar_name(self):
    return '<mixed>'

Mixed.I = Mixed()

class Sequence(Base):
  """ A sequence of values of given types """

  def __init__(self, *types):
      self.types = list(map(type_from_type, types))
      self._grammar = pp.And([i.grammar for i in self.types]).setParseAction(lambda x: tuple(x))

  def _validate(self, value):
      if not isinstance(value, tuple) or len(value) != len(self.types):
          breakpoint()
          return f'A tuple of {len(self.types)} values is required'
      for i,j in zip(self.types, value):
          out = i.validate(j)
          if out is not True: return out
      return True

  def grammar_name(self):
      return "(" + join( (f'{i}:{j.grammar_name()}' for i,j in zip(self.names, self.types) ) ) + ")"

Mixed.I = Mixed()

class Table(Base):
  """ Table with named columns, e.g.
      IQ     IREFQ       IMQ       NOQ  ITOQ  CONC
       1         1         1         1     1 1.000
       2         2         2         1     2 1.000
  """
  def __init__(self, **kwargs):
      self.names = list(kwargs.keys())
      self.sequence = Sequence( *kwargs.values() )

      grammar = pp.Suppress(pp.And(self.names) + pp.lineEnd)
      grammar += pp.OneOrMore(self.sequence.grammar() + pp.Suppress(pp.lineEnd))
      grammar.setParseAction( lambda x:  np.array(x.asList(), self.numpy_type))
      self._grammar = grammar

  @functools.cached_property
  def numpy_type(self):
      return list(zip(self.names, (i.numpy_type for i in self.sequence.types)))

  def grammar_name(self):
      data = " ".join( (f'{i}:{j.grammar_name()}' for i,j in zip(self.names, self.sequence.types) ) )
      return f"<TABLE of {data}>"


integer = Integer.I
unsigned = Unsigned.I
boolean = Bool.I
flag = Flag.I
real = Real.I
string = String.I
mixed = Mixed.I
