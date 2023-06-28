""" Grammar types for array, sequences and tables """
from collections import namedtuple
import numpy as np
import pyparsing as pp
import itertools
import copy

from ..grammar import generate_grammar, delimitedList, \
                      line_end, White
from .grammar_type import GrammarType, compare_numpy_values, type_from_type, type_from_default_value
from ..decorators import add_to_signature, cached_property
from .basic import Integer, Real, Unsigned

class Array(GrammarType):
  """ A (numpy) array of values of one type """

  delimiter=White(' \t').suppress()
  delimiter_str = ' '
  array_access = True

  def __init__(self, type, default_value=None,
               length=None, max_length=None, min_length=None,
               as_list=False, format=None,
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
    if format is not None:
       self.type = self.type.copy()
       self.type.format = format
    self.as_list = as_list
    super().__init__(default_value=default_value, **kwargs)
    self.min_length = min_length or length
    self.max_length = max_length or length
    with generate_grammar():
      grammar = self.type.grammar()
      grammar = delimitedList(grammar, self.delimiter)
      self._set_convert_action(grammar)
      grammar.setName(self.grammar_name())

    self._grammar = grammar

  def _set_convert_action(self, grammar):
    if self.as_list:
      if callable(self.as_list):
        grammar = grammar.setParseAction(lambda x: self.as_list(x.asList()))
      else:
        grammar = grammar.setParseAction(lambda x: [x.asList()])
    else:
      grammar.setParseAction(lambda x: self.convert(x.asList()))

  def __str__(self):
    if self.min_length == self.max_length:
       if self.min_length:
          length = f' of length {self.min_length}'
       else:
          length = ''
    else:
       if self.min_length is not None:
          length='{self.min_length}<=n'
       else:
          length='n'
       if self.max_length is not None:
          length+=f'<=self.max_length'
       length=' with length '
    return f"Array(of {self.type}{length})"

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

  def _validate(self, value, why='set'):
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
          self.type.validate(v, why='set')
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
       if not hasattr(value, '__iter__'):
           value = [ value ]
           ln=1
       else:
           ln=len(value)
       value = [ self.type.convert(i) for i in value ]
       out = np.asarray(value)
       return out

    return value

  is_the_same_value = staticmethod(compare_numpy_values)

class SetOf(Array):
  """ Set of values of the same type. E.g. {1,2,3} """

  delimiter = pp.Suppress(pp.Literal(',') | pp.Literal(';') | White(' \t')).setName('[,; ]')
  delimiter_str = ','

  @add_to_signature(Array.__init__)
  def __init__(self, type, *args, **kwargs):
    kwargs.setdefault('prefix', '{')
    kwargs.setdefault('postfix', '}')
    super().__init__(type, *args, **kwargs)

  def transform_grammar(self, grammar, param_name=False):
    return grammar | self.type.grammar(param_name).copy().addParseAction(lambda x: np.atleast_1d(x.asList()))

  def copy_value(self, value):
      return copy.deepcopy(value)

class Complex(SetOf):
  array_access = False

  @add_to_signature(SetOf.__init__)
  def __init__(self, *args, **kwargs):
    super().__init__(Real.I, *args, as_list=complex, length=2, **kwargs)

  def convert(self, value):
    return complex(value)

  def _validate(self, value, why='set'):
    return isinstance(value, (complex, np.complexfloating)) or 'A complex value required, {value} given.'
  def _grammar_name(self):
    return '{complex (as 2 reals)}'

  def _string(self, val):
    return real._string(val.real) + ' ' + real._string(val.imag)

  __str__ = GrammarType.__str__

class Sequence(GrammarType):
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
      if isinstance(format, dict):
         format = { type_from_type(k):v for k,v in format.items() }
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

  def _validate(self, value, why='set'):
      if not isinstance(value, (self.value_type)) or len(value) != len(self.types):
          return f'A tuple of {len(self.types)} values is required'
      for i,j in zip(self.types, value):
          out = i.validate(j, why=why)
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

  is_the_same_value = staticmethod(compare_numpy_values)

class Table(GrammarType):
  """
  Table, optionaly with named columns, e.g.

    ::text

      IQ     IREFQ       IMQ       NOQ  ITOQ  CONC
       1         1         1         1     1 1.000
       2         2         2         1     2 1.000

  """
  array_access = True
  name_in_grammar = False

  def __init__(self, columns=None,
                     header=None, free_header=False,
                     format = {float: '>22.14', None: '>16'}, format_all=True,
                     numbering=None, numbering_label=None, numbering_format=True,
                     prefix=None, postfix=None, length=None,
                     row_condition=None, flatten=False,
                     default_values=False,
                     named_result = False, **kwargs):
      if columns is None:
         columns = kwargs
         kwargs = {}
      super().__init__(prefix=None, postfix=None, **kwargs)
      if isinstance(columns, dict):
         self.names = list(columns.keys())
         columns = columns.values()
      else:
         self.names = None
      if header is None:
         header = self.names
      self.sequence = Sequence( *columns, format=format, format_all=format_all, condition = row_condition, default_values=default_values )
      self.header = header
      self.flatten = flatten
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
      if self.flatten:
          grammar.addParseAction(lambda x: x[0].ravel())
      return grammar

  def _string(self, data):
      out = []
      if self.header:
         names = ((i[1] if isinstance(i, tuple) else i) for i in self.names)
         formated = (
                 t.format_string(n) \
                 for n,t in zip(self.names, self.sequence.types)
                 )
         header = ' '.join(formated)
         if self.numbering:
            header = self.numbering.format_string(self.numbering_label or '') + ' ' + header
         out.append(header)
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

  def _validate(self, value, why='set'):
      if not isinstance(value, np.ndarray):
         return f"Numpy array as a value required {value.__class__} given"
      dtype = self.numpy_type
      dim = 1 if isinstance(dtype, list) or self.flatten else 2
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

  def number_of_columns(self):
      return len(self.sequence.types)

  def zero_data(self, length):
      """ Return array of zeros with the given number of rows and
          with the dtype of the table
      """
      dtype = self.numpy_type
      if isinstance(dtype, list):
         return np.zeros(length, dtype)
      else:
         return np.zeros((length, self.number_of_columns()), dtype)

  def grammar_name(self):
      if self.names:
        data = " ".join( (f'{i}:{j.grammar_name()}' for i,j in zip(self.names, self.sequence.types) ) )
      else:
        data = self.sequence.grammar_name()
      return f"<TABLE of {data}>"

  is_the_same_value = staticmethod(compare_numpy_values)

set_of_integers = SetOf(Integer.I)
""" A standard grammar type instance for array of integers (of any length, used by variant types) """

set_of_reals = SetOf(Real.I)
""" A standard grammar type instance for array of reals (of any length, used by variant types) """

complex_number = Complex.I = Complex()
""" A standard grammar type instance for complex numbers """
