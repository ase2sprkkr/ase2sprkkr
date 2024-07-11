""" Grammar types for array, sequences and tables """
from collections import namedtuple
import numpy as np
import pyparsing as pp
import itertools
import copy
from typing import List, Optional, Union

from ..grammar import generate_grammar, delimitedList, \
                      line_end, White
from .grammar_type import GrammarType, compare_numpy_values, type_from_type, type_from_default_value
from ..decorators import add_to_signature, cached_property
from .basic import Integer, Real, Unsigned, real


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
           raise ValueError("It is not possible for an Array to provide default_value both in 'default_value' and in 'type' argument")
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
          length+=f'<={self.max_length}'
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
       if np.__version__ >= '1.23' or self.type.numpy_type != object:

           def validate(v):
               self.type.validate(v)
               return v

           value = ( validate(self.type.convert(i)) for i in value)
           out = np.fromiter(value, dtype = self.type.numpy_type, count=ln)
       else:
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
               default_values=False, names=None, as_list=False, **kwargs):
      super().__init__(**kwargs)
      if names:
         self.names = names if isinstance(names, dict) else {name:i for i,name in enumerate(names)}
         self.value_type = namedtuple("_".join(names), names)
         self.value_constructor = self.value_type
         self.value_constructor = self.value_type
      else:
         self.names = None
         if as_list:
            self.value_type = list
            self.value_constructor = lambda *x: [[*x]]
         else:
            self.value_type = tuple
            self.value_constructor = lambda *x: [ tuple(x) ]
      if isinstance(format, dict):
         format = { type_from_type(k):v for k,v in format.items() }
      if isinstance(format, (str, dict)):
        format = itertools.repeat(format)
      self.types = [ type_from_default_value(i, dfs, format_all=format_all) for i,dfs in zip(types, format) ]
      if allowed_values and not isinstance(allowed_values, set):
         allowed_values = set(allowed_values)
      self.allowed_values=allowed_values
      self.default_values=default_values

  def __len__(self):
      return len(self.types)

  def __bool__(self):
      return True

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
          i.validate(j, why=why)
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
      return " ".join( (f'{j.grammar_name()}' for j in self.types) )

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


class SpecialColumn:
    """ A special column for :class:`Table` """
    def __init__(self, table, column, label, format, *, default_format='<4'):
        if column.__class__ is str:
           label=column
           column=True
        column = Unsigned.I if column is True else column
        if column:
           if format and not (format is True and column.format):
               if format is True:
                  format = default_format
               column = column.copy()
               column.format = format
        self.column = column
        self.label = label
        self.format = format

    def __bool__(self):
        return bool(self.column)

    def add_grammar(self, grammar):
        return self.column.grammar() + grammar

    def add_header_grammar(self, grammar):
        if self.label:
            grammar = pp.CaselessKeyword(self.label) + grammar
        return grammar

    def format_string(self, val):
        if not self.column:
            return ''
        return self.column.format_string(val)

    def header(self):
        return self.format_string(self.label or '') + ' '

    def __repr__(self):
        if not self.column:
            return "<Special column not present>"
        return f"<Special column{' ' if self.label else ''}{self.label} of type {self.column}>"


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

  def __init__(self, columns:List[GrammarType]=None,
                     header:Optional[bool]=None, free_header=False,
                     format = {float: '>22.14', None: '>16'}, format_all=True,
                     numbering:Union[str,bool,GrammarType]=None, numbering_label=None, numbering_format=True,
                     grouping=False, grouping_label=None, grouping_format=True,
                     prefix=None, postfix=None, length=None,
                     row_condition=None, flatten=False,
                     default_values=False,
                     named_result=None,
                     group_size=None, group_size_format="{:<12}{}",
                     groups_as_list=None,
                     **kwargs):
      """
      Parameters
      ----------
      columns
        List of GrammarTypes that describes the columns of the table
      header
        Whether table will have header. Default None means True, if the names are specified
      free_header
        Do not require the exact content of the header, just skip one line during parsing
      format
        Pass this argument to the :class:`Sequence` describing the line
      format_string
        Pass this argument to the :class:`Sequence` describing the line
      numbering
        There will be one extra column on the begining of the table, with numbering starting from one.
        String argument means the label of the column.
      numbering_label
        The label of the numbering column (if not given in numbering).
      numbering_format
        Format for the numering column
      grouping
        If True, the data are not one table, but list of tables.
        Numbering then numbers the tables, not the rows of the tables.
        There will be one extra column which numbers rows within the table.
      grouping_label
        Similiar as `numbering_label`
      grouping_format
        Similiard as `numbering_format`
      prefix
        Prints the prefix before the table
      postfix
        Prints the postfix after the table
      length
        The table length (or the number of groups)
      row_condition
        The condition, that each row of the table should have satisfied
      flatten
        The resulting table will have one dimension. Is meaningful for table without
        named columns, with a same type of all columns
      named_result
        The resulting table will have numpy structured dtype with named columns.
        Default None means autodetection: it behaves as True if necessary, i.e.:
        if there are at least two columns with different types.
      group_size
        If grouping, there will be another line after the header (SPRKKR format of tables),
        which contains ``{NAME}     {VALUE}``. The name should be specified in group_size.
        The value is the number of rows, which all subtables are required to have.
      group_size_format
        Format for the group_size line.
      groups_as_list
        If True - groups are contained in list
        If False - groups are contained in np.ndarray
        If None - True if group_size is defined (and thus if it is possible)
      kwargs
        Columns and their names can be assigned as kwargs, e.g.
        ``column1_name = float, column2_name = int, ...``
      """
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
         header = bool(self.names)
      self.header = header

      if group_size:
          if not grouping:
              grouping=True
      self.group_size = group_size
      self.group_size_format = group_size_format

      self.sequence = Sequence(*columns, format=format, format_all=format_all, condition = row_condition, default_values=default_values)
      self.flatten = flatten
      self.free_header = free_header

      self.named_result = named_result
      self.length = length
      self.numbering = SpecialColumn(self, numbering, numbering_label, numbering_format)
      self.grouping = SpecialColumn(self, grouping, grouping_label, grouping_format)
      self.groups_as_list = bool(self.grouping) and (
          groups_as_list if groups_as_list is not None else not group_size
      )
      if self.grouping and not self.groups_as_list and not group_size:
          raise ValueError("Groups have to be returned as list, if there is not fixed number"
              " of items in a group")

  def special_columns(self):
      if self.grouping:
          yield self.grouping
      if self.numbering:
          yield self.numbering

  def _grammar(self, param_name=False):
      line = self.sequence.grammar(param_name)
      cols = 1
      for i in self.special_columns():
          cols+=1
          line = i.add_grammar(line)
      grammar = delimitedList(line, line_end)

      if self.group_size:

          def set_g_size(x):
              grp_size.group_size = x[0]

          grp_size = Unsigned.I.grammar().copy().setParseAction(set_g_size)
          grammar = pp.Suppress(pp.CaselessKeyword(self.group_size) + grp_size + "\n") + grammar

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
                    elif n != '':
                      yield pp.CaselessKeyword(n)
             header = pp.And(list(names()))
             for i in self.special_columns():
                header = i.add_header_grammar(header)
             grammar = pp.Suppress(header + pp.lineEnd) + grammar

      def data_numbering(s, loc, x):
          numbers = x[::2]
          datas = x[1::2]
          if not numbers == [*range(1, len(numbers) + 1)]:
             raise pp.ParseException(s, loc, 'The first column should contain row numbering')
          return tabelize(datas)

      def data_grouping(s, loc, data):

           out = []
           curr = []
           c = 0
           if len(data) % (2 * grp_size.group_size):
              raise pp.ParseException(s, loc, "The data in the table should be grouped "
                                           f"by groups of length {grp_size.group_size}.")
           for i in range(0, len(data), 2):
               c+=1
               if c>grp_size.group_size:
                   out.append(tabelize(curr))
                   curr = []
                   c = 1
               if data[i] != c:
                  raise pp.ParseException(s, loc, "The value in the first column should "
                        "be a row numbering within a group of a length {grp_size.group_size}."
                       f"{data[i]} encountered, {c} expected.")
               curr.append(data[i + 1])
           out.append(tabelize(curr))
           if self.groups_as_list:
                return [ out ]
           else:
                return np.asarray(out, dtype = self._numpy_type)

      def data_numbering_grouping(s, loc, data):
           if len(data) == 0:
               return data

           def numbering_error(i, expected):
               raise pp.ParseException(s, loc, "The value in the first column should "
                     "be a monotone numbering of a groups. "
                     "On the {i//3+1}th row, {data[i]} encountered, {expected} expected.")
           val = data[0]
           if val != 1:
                numbering_error(0, "1")
           out = []
           curr = []
           grp = 0
           for i in range(0, len(data), 3):
               if data[i] !=val:
                  if data[i] != val + 1:
                       numbering_error(i, f"'{val}' or '{val+1}'")
                  out.append(tabelize(curr))
                  curr = []
                  val = data[i]
                  grp = 0
               grp += 1
               if data[i + 1] != grp:
                  raise pp.ParseException(s, loc, "The value in the second column should "
                        "be a row numbering within a group. "
                       f"{data[1+1]} encountered, {grp} expected.")
               curr.append(data[i + 2])
           out.append(tabelize(curr))
           if self.group_size:
              for i,g in enumerate(out):
                  if len(g) != grp_size.group_size:
                      raise pp.ParseException(s, loc, f"The group {i//3 + 1} should have "
                              f"{grp_size.group_size} members, it has {len(g)}.")
           if self.groups_as_list:
                return [ out ]
           else:
                return np.asarray(out, dtype = self._numpy_type)

      def tabelize(x):
          out = np.array(x, dtype=self._numpy_type)
          if self.flatten:
              out = out.ravel()
          return out

      if self.numbering:
          if self.grouping:
              grammar.addParseAction(data_numbering_grouping)
          else:
              grammar.addParseAction(data_numbering)
      elif self.grouping:
          grammar.addParseAction(data_grouping)
      else:
          grammar.addParseAction(tabelize)

      return grammar

  def _string(self, data):
      out = []
      if self.header:
         names = ((i[1] if isinstance(i, tuple) else i) for i in self.names)
         formated = (
                 t.format_string(n)
                 for n,t in zip(names, self.sequence.types)
                    )
         header = ' '.join(formated)
         for i in self.special_columns():
             header = i.header() + header
         out.append(header)

      if self.group_size:
         out.append(self.group_size_format.format(self.group_size, len(data[0]) if data is not None and len(data) else 1))

      def line(row, prefix):
          line = self.sequence.string(row)
          if prefix:
              line = prefix + line
          out.append(line)

      def unflatten(val):
          if val.ndim == 1 and not isinstance(self._numpy_type, list):
              ln = len(self.sequence)
              val = val.reshape(len(val) // ln, ln)
          return val

      if self.grouping:
          for n, group in enumerate(data, 1):
              for g, d in enumerate(unflatten(group), 1):
                  line(d, self.grouping.format_string(n) + self.numbering.format_string(g))
      else:
          for n, d in enumerate(unflatten(data),1):
              line(d, self.numbering.format_string(n))

      return '\n'.join(out)

  def _validate(self, value, why='set'):
      def validate(value):
          if not isinstance(value, np.ndarray):
             return f"Numpy array as a value required {value.__class__} given"
          dtype = self._numpy_type
          dim = 1 if isinstance(dtype, list) or self.flatten else 2
          if len(value.shape) != dim:
             return f"The array should have dimension={dim}, it has dimension {len(value.shape)}"
          if value.dtype != dtype:
             return f"The data type of the value should be {dtype}, it is {value.dtype}"
          if dim==2 and value.shape[1] != len(self.sequence.types):
             return f"The array is required to have {len(self.sequence.types)} columns, it has {value.shape[1]}"
          return True

      if not self.grouping:
          out = validate(value)
          if out is not True:
              return out
      else:
          if not isinstance(value, (list, np.ndarray)):
              return "The grouped table accepts list of table data as a value, "\
                     f"{value.__class__} have been given."
          for i,d in enumerate(value):
              out = validate(d)
              if out is not True:
                  return f"The {i+1}th table has invalid data: {out}"
          size = None
          if self.group_size:
              for i, d in enumerate(value):
                  if size is None:
                      size = len(d)
                  else:
                      if size != len(d):
                         return f"Size of the {i+1}th table should be same as of the first one, i.e. {size}"

      if self.length is not None and self.length != len(value):
          return f"The array is required to have {self.length} rows, it has {len(value)}"
      return True

  def convert(self, value):
      if self.groups_as_list:
          value = [ np.asarray(v, dtype = self._numpy_type) for v in value ]
      return np.asarray(value, dtype = self._numpy_type)

  @cached_property
  def numpy_type(self):
      if self.groups_as_list:
          return object
      return self._numpy_type

  @cached_property
  def _numpy_type(self):
      types = self.sequence.types
      if not self.named_result:
         dtype = types[0].numpy_type
         for t in types[1:]:
             if t.numpy_type != dtype:
                 if self.named_result is False:
                      raise ValueError("The table should not have structured dtype, but has columns have different types")
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

complex_number = Complex.I = Complex()        # NOQA E741
""" A standard grammar type instance for complex numbers """
