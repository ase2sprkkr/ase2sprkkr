""" Basic types for GrammarTypes and usefull functions """
import copy
from typing import Union, Any, Callable, Optional, Type, Dict, List
import functools
from collections.abc import Hashable
import numpy as np
import pyparsing as pp
import inspect
from .. import grammar_types

from ..decorators import cached_class_property, cache, \
                         add_called_class_as_argument, cached_property
from ..alternative_types import normalize_type, allowed_types
from ..grammar import generate_grammar
from ..formats import full_format_for_string


class GrammarType:
  """ Base class for definition of configuration option types

      A type without value (e.g. Separator) are just syntactical
      elements in the potentials file, that do not carry an information.
      Such elements do not yields (name, value) pair during parsing the file.

      Do not confuse this with GrammarType.missing_value functionality.
      Missing_value is just the opposite: missing_value can be ommited in the file
      (or even the absence of the name in the file carry the information, that
      the Flag is False), but the name-value tuple of such Type is present
      in the parse result. On the other hand, has_value = False is in the file, but
      not in the result.

      **The functions called during...**

      ::

        User input:  convert, validate

        Output: string -> _string

        Parsing: parse -> ( <_grammar parse actions>, validate(why='parse') )

  """

  has_value = True

  name_in_grammar = True
  """ Default value for ValueDefinition.name_in_grammar.
      Some types (e.g. Tables) commonly have no name (are identified
      by its position in the potential file) -- such type could redefine
      this class property."""

  default_value = None
  """ Default value for the given type. It can be overriden for particular instances
  in the constructor (or just by setting the attribute of an instantiated object). """

  numpy_type = object
  """ The numpy dtype of the array, that contains values of this type (see e.g. :class:`Array`).
      The default type ``object`` can and should be redefined in the descendatns. """

  array_access = False
  """ The value of this type can be accessed as array """

  is_independent_on_the_predecessor = False
  """ Options of most grammar types do not identify themselves, so they have to be either identified by their names, or if name is not given, by their predecessors. Hoewever, e.g. keyword arguments can be identified just by their value."""

  def __init__(self, prefix:Union[str,None]=None, postfix:Union[str,None]=None,
                     format:str='', after_format:Optional[str]=None,
                     default_value:Any=None,
                     condition:Union[Callable[[Any], Union[bool,str]],None]=None,
                     after_convert:Union[Callable[[Any], Any],None]=None,
                     description=''):
      """
      Create the object.

      Parameters
      ----------
      prefix
        The string, that will be printed before the value

      postfix
        The string, that will be printed after the value

      format
        The (python) format string, that will be used for outputing the value.
        The format is passed as format argument to ``str.format`` routine.

      after_format
        In some cases, the additional formating is required after converting to the string
        and adding postfix/prefix.

      default_value
        The default value of the options of this type. ``None`` means no default value.

      condition
        Function, that check the validity of the value. It should return ``True`` for a valid
        value, and ``False`` or string for invalid. The string is interpreted as an error message
        that explains the invalidity of the value.

      after_convert
        Function, that - if it is given - is applied to the (entered or parsed) value. The function
        is applied on the result of the
        :meth:`convert<ase2sprkkr.common.grammar_types.GrammarType.convert>` method
      """

      self.prefix = prefix
      """ The string, that will be printed before the value """
      self.postfix = postfix
      """ The string, that will be printed after the value """
      self._format = format
      self.after_format = after_format if not after_format or '{' in after_format else \
                         f'{{:{after_format}}}'
      """ The (python) format string, that will be used for printing the value.
        The format is passed as format argument to ``str.format`` routine.  """
      self.condition = condition
      if after_convert is not None:
         self.convert = lambda v: \
              after_convert(self, self.__class__.convert(self, v))

      """ Some subclasses has default_value defined via read-only property. """
      if default_value is not None:
         self.default_value = self.convert(default_value)
      self._description = description

  def __str__(self):
      return self.__class__.__name__

  @cached_property
  def format(self):
      """ Return the resulting format string, applying the prefix and postfix """
      if not self._format and not self.prefix and not self.postfix:
          return None
      out = self._format or ''
      if '{' not in out:
            out=f'{{:{out}}}'
      escape = lambda x: x.replace('{','{{').replace('}', '}}')
      if self.prefix:
        out=escape(self.prefix) + out
      if self.postfix:
        out+=escape(self.postfix)
      return out

  @format.setter
  def format(self, v):
      self._format = v
      if 'format' in self.__dict__:
          del self.__dict__['format']

  @staticmethod
  def is_the_same_value(a,b):
    """ Comparison function for the values of "this type".

    Not all values (e.g. numpy arrays) can be compared by equal sign,
    so this function has to be used for comparison of the values.
    """
    return a == b

  @cache
  def grammar(self, param_name:str=False):
    """ Return a pyparsing grammar for the type

    Parameters
    ----------

    param_name
      The name of the value, that can be assigned to the generated grammar element.
    """
    grammar = self._grammar

    if isinstance(self._grammar, pp.ParserElement):
       grammar = pp.Forward()
       grammar << self._grammar
    else:
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
             self.validate(x[0], why='parse', param_name=param_name)
           except ValueError as e:
             raise pp.ParseException(s, loc, str(e) + '\nValidating of the parsed value failed') from e
           return x

       grammar.addParseAction(validate)
    grammar.grammar_type = self
    return grammar

  def parse(self, str, whole_string=True):
    """
    Parse the string, return the obtained value.
    """
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
    if not isinstance(self.grammar, pp.ParserElement):
       return self.__class__.__name__
    return str(self.grammar)

  def transform_grammar(self, grammar, param_name=False):
    """ The chance for the resulting class to alter the resulting prefixed grammar """
    return grammar

  def missing_value(self):
    """ Is the configuraion value a flag? I.e., can be =<value> ommited
    in the configuration?

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

  def validate(self, value, param_name='<Unknown>', why:str='set'):
    """ Validate either the pyparsing result or a user given value.

    Do not override this method in subclasses for the validation implementation,
    this method calls :meth:`_validate`, which should contain the actual validation

    Parameters
    ---------
    value : mixed
      Value to be validated.

    param_name : str or callable
      Parameter name to be used in possible throwed exception (Optional).
      If it is callable, it should be a function that returns the param_name.

    why
      Possible values are:

      ``set``
         validation value setted by user (in rare cases, such value can be incomplete
         and requires `completing` during ``set_from_atoms`` call before saving the output)
      ``parse``
         validation during parsing input file, checks enforced
         by the grammar can be skipped
      ``save``
         validation before saving the values
    """
    try:
      err = self._validate(value, why)
    except ValueError as err:
      self._valueError(value, err, param_name)
    if err is not True:
      self._valueError(value, err, param_name)
    if self.condition:
      err = self.condition(value)
      if err is not True:
        self._valueError(value, err, param_name)
    return True

  def _validate(self, value, why='set'):
    """ Return error message if the value is not valid. """
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

  def convert(self, value):
    """ Convert a value from user to the "cannonical form" """
    return value

  def _string(self, val):
    """
    Convert the value to the ouput.

    The :meth:`string` apply format and do some additional transformation (add prefix, postfix etc.),
    so the actual way how to convert the value for the output should be here. """
    return val

  def string(self, val):
    """ Convert the value to the string according to the class definition.

    Before redefining this method, you should consider, whether :meth:`_string` method could be
    redefined instead. Otherwise, you should call :meth:`apply_format` in the redefined method.
    to retain the common functionality (as adding prefix or postfix to the resulting
    string).
    """
    val = self._string(val)
    val = self.apply_format(val)
    return val

  def apply_format(self, val):
    """ Apply format to the outputed value. """
    if self.format:
       val = self.format.format(val)
    else:
       val = str(val)
    if self.after_format:
      val = self.after_format.format(val)
    return val

  def format_string(self, val):
    """ Format string in a similiar manner as a value.
        It is usefull for simple types, where header of a table
        should be formatted in the same way.
        For complex types it may not give a reasonable results.
    """
    if self.format:
       out=full_format_for_string(self.format).format(val)
    else:
       out=str(val)
    if self.after_format:
       out = self.after_format.format(out)
    return out

  def write(self, f, val):
    """ Output the value to the stream (in the propper format). """
    f.write(self.string(val))

  def print(self, val):
    """ Output the value to stdout (in the propper format). """
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

  def additional_description(self, prefix='') -> str:
    """ If the description of the type does not fit on one line,
    this method should return

    Returns
    -------
    additional_description
      The additional description (e.g. possible choices) of the type. Multiline string.
    """
    out = self._description
    if prefix and out:
       out = out.replace('\n', '\n' + prefix)
    return out

  def type_validation(self, value, types:Union[List[Type], Type], typename:Optional[str]=None):
    """
    Parameters
    ----------
    value: mixed
      Value to be checked

    types
      The required type or types. If more types is given, it is sufficient, if the value is of
      any of given types.

    Returns
    -------
    error_message: Union[str, bool]
      The function returns either False, if the value is ok, or string containing an error
      message describing the error.

    """
    if isinstance(value, types):
        return True
    if not typename:
       typename = types
    typename=str(typename)
    n = 'n' if typename[0] in ['a','e','i','o','u'] else ''
    return f"A{n} <{typename}> value is required, a value {value} of type {value.__class__} have been given"

  def copy_value(self, value):
    return value

  def used_in_definition(self, definition):
      pass

  def added_to_container(self, definition):
      pass


@add_called_class_as_argument
def add_to_parent_validation(validation):

    @functools.wraps(validation)
    def wrapped(cls, self, value, why='set'):
        out = super(cls, self)._validate(value, why)
        if out is not True:
           return out
        return validation(self, value, why)

    return wrapped


class TypedGrammarType(GrammarType):

  @cached_class_property
  def datatype(cls):
      """ The (primary) type of the value. Redefine it in the descendants, if it is needed. """
      return cls.numpy_type

  @cached_class_property
  def allowed_types(cls):
      """ Allowed alternative types, that will be converted to the 'primary' datatype. """
      dt = cls.datatype
      return allowed_types.get(dt, (dt, ))

  def convert(self, value):
      if isinstance(value, self.datatype):
         return value
      for i in self.allowed_types:
          if isinstance(value, i):
             return self.datatype(value)
      return value

  @cached_class_property
  def datatype_name(cls):
      return cls.__name__.lower()

  def _validate(self, value, why='set'):
      return self.type_validation(value, self.allowed_types, self.datatype_name)


def type_from_type(type, format:Union[str,Dict]='', format_all:bool=False, type_map:Dict={}):
  """ Guess and return the grammar element (GrammarType class descendatnt) from a python type. E.g. int => Integer.

      The given format can be optionally set to the returned grammar element.

      Parameters
      ----------
      type: A python type or GrammarType
        A type to be converted to a grammar type (GrammarType class descendant)

      format
        The format to be applied to the resulting class. If dict is given, see :func:`format_for_type`
        for the way how the format is determined

      format_all
        If False (default), the format is not applied, if instance of GrammarType is given as
        the type parameter. Otherwise, a copy of the input type with the applied format is returned

      type_map
  """
  if isinstance(type, dict):
      return grammar_types.Keyword(type)

  type_from_type_map = grammar_types.type_from_type_map

  if isinstance(type, GrammarType):
     if format_all:
        type = type.copy()
        type.format = format_for_type(format, normalize_type(type.numpy_type))
     return type

  if isinstance(type, Hashable):
    type = normalize_type(type)
    if type in type_map:
       type = type_map[type]
    elif type in type_from_type_map:
       type = type_from_type_map[type]
    else:
       return type

    format = format_for_type(format, type)
    if format:
        type = type.copy()
        type.format = format
  return type


def type_from_value(value, type_map={}):
  """ Gues the grammar type from a python value.

  ..doctest::
  >>> type_from_value(2)
  <Integer>

  >>> type_from_value(2.0)
  <Real>
  """

  type_from_set_map = grammar_types.type_from_set_map

  if isinstance(value, recognized_set_types):
     return type_from_set_map[normalize_type(value[0].__class__)] if len(value) else grammar_types.Integer.I
  if isinstance(value, str):
     try:
        grammar_types.String._grammar.parseString(value, True)
        return grammar_types.String.I
     except Exception:
        return grammar_types.QString.I
  if isinstance(value, dict):
      return grammar_types.Keyword(value)
  type = type_from_type(value.__class__, type_map=type_map)
  if type is value.__class__:
     raise ValueError(f'Cannot determine grammar type from value {value}')
  return type.__class__(default_value = value)


def type_from_default_value(value, format='', format_all=False, type_map={}):
   """ Guess the grammar type from a value, that will become the default value of the grammar type.

   It has to create a new object instance, as it has to set the default
   value property of the returned object. An (output) format can be applied to the
   resulting grammar type

   Grammar types passed as types are left as is, unless format_all flag is set.
   """
   if inspect.isclass(value) or isinstance(value, GrammarType):
      return type_from_type(value, format=format, format_all=format_all, type_map={})
   ptype = normalize_type(value.__class__)
   gtype = type_from_type(value.__class__, type_map=type_map).__class__
   return gtype(default_value = value, format=format_for_type(format, ptype))


def compare_numpy_values(a,b):
    """ The numpy arrays cannot be compared by =, that's why this method.
    However, the method is still far from to be perfect, it can not
    compare nested numpy arrays.
    """
    return np.array_equal(a,b)


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


recognized_set_types = ( list, tuple, np.ndarray )
""" The types, that are recognized as 'list of values' and so that will
be accepted as values for array_like type (e.g. :class:`Array` or :class:`SetOf`). """
