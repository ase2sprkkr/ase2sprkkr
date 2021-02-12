import pyparsing
from pyparsing import Word, Suppress, delimitedList
ppc = pyparsing.pyparsing_common

""" s of Configuration Values """

class Base:
  """ Base type for configuration types """

  def __str__(self):
    return self.__class__.__name__

  def grammar(self):
    """ Return a pyparsing grammar for the type """
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
      self.__valueError(value, err)
    return True

  def _validate(self, value):
    """ Return error message if the value is not valid """
    return True

  def _valueError(self, value, error_message=False):
    if error_message is False:
       error_message = 'invalid value'
    raise ValueError("Value {} for parameter {} is not valid: {}".format(value, self.param_name, message))

  def read(self, token, parameterName='<Unknown>'):
    """ Transform pyparsing token to a validated value """
    val = self.convert_token(token)
    self.validate(val)
    return val

  def convert_token(self, value):
    """ Convert token returned from pyparsing """
    return value[0]

  def convert(value):
    """ Convert a value from user to a "cannonical form" """
    return value

  def write(f,val):
    f.write(val)

class Unsigned(Base):

  _grammar = ppc.integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(value, param_name):
    if not isinstance(value, int): return "Integer value required"
    return value >= 0 or "Positive value required"

  def grammar_name(self):
    return '+-<int>'

class Integer(Base):

  _grammar = ppc.signed_integer.copy().setParseAction(lambda x:int(x[0]))

  def _validate(self, value):
    return isinstance(value, int) or "Integer value required"

  def grammar_name(self):
    return '<int>'


class Real(Base):

  _grammar = ppc.fnumber.setParseAction(lambda x: float(x[0]))


  def _validate(self, value):
    return isinstance(value, float) or "Float value required"

  def grammar_name(self):
    return '<float>'

class String(Base):

  _grammar = Word(pyparsing.printables).setParseAction(lambda x:x[0])

  def _validate(self, value):
    return isinstance(value, str) or "String value required"

  def grammar_name(self):
    return '<str>'

class Keyword(Base):

  def __init__(self, *keywords):
    self.keywords = kewords

  def _validate(self, value):
    return value in keywords or "Required one of [" + "|".join(self.keywords) + "]"

  def grammar(self):
    return pp.MatchFirst((pp.CaselessKeyword(i) for i in self.keywords)).setParseAction(lambda x: x[0].upper())

class Flag(Keyword):

  def __init__(self, flag):
     super().__init__(flag)

  def _validate(self, value):
      return value is True or "This is Flag with no value, please set to True to be present or "
    

class SetOf(Base):
  """ Set of values, e.g. {1,2,3} """
  def __init__(self, type, length=None, max_length=None):
    self.min_length = length
    self.max_length = max_length or length
    self.type = map_type(type)

  def __str__(self):
    return "SetOf({})".format(str(self.type))

  def grammar_name(self):
      return '{' + self.type.grammar_name() + ",...}"

  def _validate(self, value):
    for i,v in enumerate(value):
        out = self.type._validate(v)
        if out is not True:
           return "Value {} in set is incorect: {}".format(i, out)
        return True

  def convert(self, value):
    if self.length is None and not isinstance(list(value)):
      return [ value ]

  def grammar(self):
    return Suppress("{") + delimitedList(self.type.grammar()) + Suppress("}")

  def convert_token(self, token):
    return list(token)

""" Map python native types to configuration value types """

map_type_map = {
    int  : Integer,
    str  : String,
    float: Real
}

def map_type(type):
  if type in map_type_map:
    return map_type_map[type]()
  return type
