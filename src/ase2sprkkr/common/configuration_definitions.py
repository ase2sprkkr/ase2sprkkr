"""
Configuration definitions are classes, that desribes
the syntax of a configuration file, or its parts
(sections or configuration options)

They are able both to parse a file, which results in an
instance of (an instance of :py:class:`ase2sprkkr.common.Configuration`,
e.g. an :py:class:`Option<ase2sprkkr.common.options.Option>` or
:py:class:`Section<ase2sprkkr.common.configuration_containers.Section>`
), or write such object to a file.
"""

from ..common.grammar_types  import type_from_type, type_from_value, GrammarType
from ..common.grammar  import delimitedList, end_of_file, generate_grammar
from ..common.misc import OrderedDict
from .configuration_containers import Section
from .options import Option, DangerousValue
from .decorators import cache

import numpy as np
import pyparsing as pp
import inspect
import itertools
import builtins
from typing import Dict, Union

#:This serves just for dealing with various pyparsing versions
_parse_all_name = 'parse_all' if \
  'parse_all' in inspect.getfullargspec(pp.Or.parseString).args \
  else 'parseAll'

def unique_dict(values):
    """ Create a dictionary from the arguments. However, raise an exception,
    if there is any duplicit key.
    Moreover, if there is key of type (a,b), it will be transformed to subdictionary.

    >>> unique_dict( [ ('x', 'y'), (('a','b'), 1 ), (('a', 'c'), 2) ] )
    {'x': 'y', 'a': {'b': 1, 'c': 2}}
    """
    out = {}
    duplicates = []

    for k, v in values:
        if isinstance(k, tuple):
           if not k[0] in out:
              o = out[k[0]] = {}
           elif not isinstance(out[k[0]], dict):
              duplicates.append(k[0])
              continue
           else:
              o = out[k[0]]
           if k[1] in o:
              duplicates.append(''.join(k))
           else:
              o[k[1]] = v
        else:
           if k in out:
              duplicates.append(k)
           else:
              out[k] = v

    if duplicates:
       duplicates = ','.join(duplicates)
       raise pp.ParseException(f"There are non-unique keys: {duplicates}")
    return out

class BaseDefinition:
   """ A base class for a configuration definition, either of an option, or of a container.
   The definition determine the type of value(s) in the configuration/problem-definition file,
   the format of the values, the way how it/they is/are read etc...
   """

   result_class = None
   """ Parsing of the data results in an instance of this class. To be redefined in descendants. """

   name_in_grammar = True
   """ Is the name of the value/section present in the file?

   When False, the option/section name is not written to the file at all. The value(s)
   of the option/section is/are then expected to be located in the file
   just after the previous one.

   By default, all options and sections… are named (in the configuration file). However,
   the attribute can be redefined in instantiated objects and/or descendant classes to change
   the behavior.
   """

   def __init__(self, name, alternative_names=None,
                is_optional=False, is_hidden=False, is_expert=False,
                name_in_grammar=None, info=None, description=None,
                write_alternative_name:bool=False,
                write_condition=None,
                result_class=None
                ):
       """
       Parameters
       ----------
        name: str
          Name of the value/section

        alternative_names: str or [str]
          Alternative names that can denotes the value

        is_optional: boolean
          If True, this section/value can be missing in the .pot/task file

        is_hidden: boolean
          Hidden values are not offered to a user, usually they are
          set by another object (and so a direct setting of their values
          has no sense)

        is_expert: boolean
          Expert values/sections are not required and they are somewhat hidden
          from the user

        name_in_grammar: boolean or None
          If False, there the name of the variable is not printed in the
          configuration file. The variable is recognized by its position.
          If None, the default class value is used

        info: str
          A short help message for the value/section. It will be the perex for description.

        description: str
           The additional informations for the users.

        write_alternative_name
           Wheter use the name or the (first) alternative name in the output.

        write_condition
           If defined, write the value, only if write_condition(the option) is True.

        result_class
           Redefine the class that holds data for this option/section.
       """
       self.name = name
       """ The name of the option/section """
       if isinstance(alternative_names, str):
          alternative_names = [ alternative_names ]
       self.alternative_names = alternative_names
       """ Alternative names of the option/section. The option/section can
       be "denoted" in the configuration file by either by its name or any
       of the alternative names.
       """
       self.is_optional = is_optional
       self.is_expert = is_expert
       self.is_hidden = is_hidden
       """ Is it required part of configuration (or can it be ommited)? """
       self.write_alternative_name = write_alternative_name
       self.write_condition = write_condition
       self.name_in_grammar = self.__class__.name_in_grammar \
                               if name_in_grammar is None else name_in_grammar
       self._info = info
       """ A short help text describing the content for the users. """
       self._description = description
       """ A longer help text describing the content for the users. """
       if result_class:
           self.result_class = result_class

   def info(self, generic:bool=True) -> str:
       """ Return short help string.

       Parameters
       ----------
       generic
         If the definition has no help specified and generic is True, return a (not saying much) generic help string
       """
       out = self._info
       if not out and generic:
          out = self._generic_info()
       return out


   def create_object(self, container=None):
       """ Creates Section/Option/.... object (whose properties I define) """
       return self.result_class(self, container)

   def grammar(self, allow_dangerous:bool=False):
       """ Generate grammar with the correct settings of pyparsing global state

       Parameters
       ----------
       allow_dangerous
        Allow dangerous values - i.e. values that do not fulfill the requirements
        for the given-option value (i.e. a type requirement or other constraints).
       """
       with generate_grammar():
         return self._grammar(allow_dangerous)

   _description_indentation = '    '
   """ Nested levels of description will be indented using this 'prefix' """

   def description(self, verbose:bool=True, show_hidden=False, prefix:str=''):
       """
       Parameters
       ----------
       verbose: bool
         If true, add also detailed documentation of all included items (e.g. members of a container)

       show_hidden: bool
         If ``False`` (default), do not show descriptions of hidden members.

       prefix
         The string, with with each line will begin (commonly the spaces for the indentation).
       """
       out = [
          prefix + self.info().replace('\n', '\n' + prefix), prefix
       ]

       #out.append(f"{prefix}Data description\n"
       #           f"{prefix}----------------")
       out.append(self.data_description(verbose, show_hidden, prefix))
       if self._description:
          out.append('\n')
          out.append(prefix + self._description.replace('\n', '\n' + prefix))
       return '\n'.join(out)

   def _tuple_with_my_name(self, expr,
                           delimiter=None,
                           has_value:bool=True,
                           is_numbered_array:bool=False):
        """ Create the grammar returning tuple (self.name, <result of the expr>)

            Parameters
            ----------
            expr
              Pyparsing expresion for the value of the option/section
            delimiter
              Pyparsing expression for the name-value delimiter
            has_value
              If False, do not add the parsed value to the results.
              This can be used e.g. for separators (see :class:`ase2sprkkr.common.grammar_types.Separator`) etc.
            is_numbered_array
              If True, the resulting grammar is in the form
              NAME[index]=....
        """
        if self.name_in_grammar:
            keyword = pp.CaselessLiteral if is_numbered_array else pp.CaselessKeyword
            name = keyword(self.name)
            if self.do_not_skip_whitespaces_before_name:
               name.leaveWhitespace()

            if self.alternative_names:
               alt_names = ( keyword(n) for n in self.alternative_names )
               if self.do_not_skip_whitespaces_before_name:
                  alt_names = ( i.leaveWhitespace() for i in alt_names )
               name = name ^ pp.Or(alt_names)
               if self.do_not_skip_whitespaces_before_name:
                  name.leaveWhitespace()
            name.setParseAction(lambda x: self.name)
            if is_numbered_array:
               name += pp.Optional(pp.Word(pp.nums), default='def')
               name.setParseAction(lambda x: (x[0], 'def' if x[1]=='def' else int(x[1])) )
            if delimiter:
              name += delimiter
            out = name - expr
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)
            out = name + expr
        if has_value:
            return out.setParseAction(lambda x: tuple(x))
        else:
            return out.suppress()

   do_not_skip_whitespaces_before_name = False

   def _get_copy_args(self)->Dict[str, str]:
       """
       Compute the dictionary that defines the attributes to create a copy of this object.

       Returns
       -------
       copy: Dict
          The returning dictionary has this structure:
          { name of the argument of the __init__ function : name of the object attribute }
       """

       if not '_copy_args' in self.__class__.__dict__:
          args = inspect.getfullargspec(self.__class__.__init__).args[1:]
          self.__class__._copy_args = { v: '_'+v if '_'+v in self.__dict__ else v for v in args if v not in self._copy_excluded_args }
       return self.__class__._copy_args

   _copy_excluded_args = ['expert']

   def can_be_repeated(self):
       """ If True, the item can be repeated in the parsed file.
       This behavior have currently the values with is_numbered_array property = True.
       This function is to be redefined in descendants
       """
       return False

   is_generated = False
   """ Generated values are computed on the fly from the other data """

class ValueDefinition(BaseDefinition):

  result_class = Option

  name_in_grammar = None

  def __init__(self, name, type=None, default_value=None, alternative_names=None,
               fixed_value=None, required=None, info=None, description=None,
               is_hidden=False, is_optional=None, is_expert=False, is_numbered_array:bool=False,
               is_always_added:bool=None,
               name_in_grammar=None, name_format=None, expert=None,
               write_alternative_name:bool=False, write_condition=None,
               result_class=None,
               ):
    """
    Definition of a configuration value.

    Parameters
    ----------
    name: str
      Name of the configuration value

    type: Optional[GrammarType|mixed]
      Configuration value data type.
      If it is set to anyting what is not derived from GrammarType, the given value is used as the default value
      and the data type is derived from it.
      If it is None, the default value have to be set using ``expert`` parameter.

    default: mixed
      Default value for the configuration option.
      Can accept callable (with option instance as an argument) - then the value will be determined
      at 'runtime' (possibly according to the other values of the section)

    alternative_names: str or [str]
      Value can have an alternative name (that alternativelly denotes the value)

    fixed_value: mixed
      If it is given, this option have a fixed_value value (provided by this parameter),
      that can not be changed by an user.
      #TODO - currently, callback (as in default_value) is not supported

    required: bool
      Required option can not be set to None (however, a required one
      can be still be optional, if it has a default values).
      If required = None, it is set to True if both the conditions are met:

       * the value is not expert
       * the optional is not True and the option has not a default_value

    is_optional: bool or None
      If True, the value can be omited, if fixed order in the section is required
      None means True just if required is False (or it is determined to be False),
      see the ``required`` parameter.

    is_hidden: bool
      The value is hidden from the user (no container.name access to the value).

    is_expert: Union[bool,mixed]
      Expert values are somewhat hidden (e.g. listed at end) from the user.
      Expert values are not exported to the result, if they are set to the
      default value.

    is_numbered_array
      Such values can contains (sparse) arrays. In the resulting ouput, the
      members of the array are in the form NAME1=..., NAME2=..., ... The default
      value for missing number can appear in the form NAME=...

    is_always_added
      If False, add the value, only if its value is not the default value.
      Default None means False for expert values, True for the others.

    name_in_grammar: bool or None
      The value in the conf file is prefixed by <name><name_value_delimiter>
      If None, the default type value (type.name_in_grammar) is used

    name_format: str or None
      The way how the name is written

    expert: Optional[mixed]
      If not None, set ``is_expert`` to True, ``default_value`` to the given value and
      ``required`` to False. Note, that also ``type`` can be determined from such given
      ``default_value``.

    write_alternative_name
       Wheter use the name or the (first) alternative name in the output.

    write_condition
       If defined, write the value, only if write_condition(the option) is True.

    result_class
       Redefine the class that holds data for this option/section
    """
    if expert is not None:
       if type is None:
          type = expert
       else:
          default_value=expert
       is_expert = True
       required = False
    elif type is None:
       raise TypeError("The data-type of the configuration value is required.")

    if is_always_added is None:
       self.is_always_added = not is_expert
    else:
       self.is_always_added = is_always_added

    if fixed_value is None:
       self.is_fixed = False
    else:
       default_value = fixed_value
       self.is_fixed = True

    if default_value is None and not isinstance(type, (GrammarType, builtins.type)):
       self.type = type_from_value(type, type_map = self.type_from_type_map)
       self.default_value = self.type.convert(type)
    else:
       self.type = type_from_type(type, type_map = self.type_from_type_map)
       self.default_value = self.type.convert(default_value) if default_value is not None else None
    assert isinstance(self.type, GrammarType), "grammar_type (sprkkr.common.grammar_types.GrammarType descendat) required as a value type"

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    if required is None:
       required = not is_expert and (not is_optional and default_value is None)

    if is_optional is None:
       is_optional = required is False

    configuration_type_name = 'OPTION'

    super().__init__(
         name = name,
         alternative_names = alternative_names,
         is_optional = is_optional,
         is_hidden = is_hidden,
         is_expert = is_expert,
         name_in_grammar = name_in_grammar,
         info=info,
         description = description,
         write_alternative_name = write_alternative_name,
         write_condition = write_condition,
         result_class = result_class
    )

    if self.name_in_grammar is None:
        self.name_in_grammar = self.type.name_in_grammar

    self.is_numbered_array = is_numbered_array
    if is_numbered_array and not self.name_in_grammar:
       raise ValueException('Numbered_array value type has to have its name in the grammar')

    self.required = self.default_value is not None if required is None else required
    self.name_format = name_format

  type_from_type_map = {}
  """ Redefine this in descendants, if you need to create different types that the defaults to be
  'guessed' from the default values """

  def enrich(self, option):
      """ The Option can be enriched by the definition, e.g. the docsting can be extended. """
      self.type.enrich(option)

  @property
  def formated_name(self):
    name = next(iter(self.alternative_names)) if self.write_alternative_name else self.name
    if self.name_format:
       return "{:{}}".format(name, self.name_format)
    return name

  def data_description(self, verbose:Union[bool,str]=False, show_hidden=False, prefix:str=''):
    """
    Return the description of the contained data type and their type.

    Parameters
    ----------
    verbose
      If ``False``, return only one-line string with a basic info.
      If ``True``, return more detailed informations.
      'verbose' means here the same thing as True

    show_hidden
      If `False``, do not show hidden members... which has no meaning for Values.

    prefix
      The string, with with each line will begin (commonly the spaces for the indentation).
    """
    out = f"{prefix}{self.name} : {self.type}"
    if callable(self.default_value):
        value = getattr(self.default_value,'__doc__', '<function>')
    else:
        value = self.get_value()
    if value is not None:
       out+=f" ≝ {value}"

    flags = []
    if self.is_optional:
       flags.append('optional')
    if self.is_hidden:
       flags.append('hidden')
    if self.is_expert:
       flags.append('expert')
    if self.is_expert == self.is_always_added:
       flags.append('always add' if self.is_expert else 'add non-default')
    if self.is_numbered_array:
       flags.append('array')
    if self.is_fixed:
       flags.append('read_only')
    if flags:
       flags = ', '.join(flags)
       out+= f"  ({flags})"

    if verbose:
       add = self.additional_data_description(prefix=prefix + self._description_indentation)
       if add:
          out+= f'\n'
          out+= add
    return out

  def additional_data_description(self, verbose=False, show_hidden=False, prefix:str='') -> str:
    """ Return the additional runtime-documentation for the configuration value.
        E.g. return the possible choices for the value, etc...

        Parameters
        ----------
        verbose
          This parameter has no effect here. See :meth:`BaseDefinition.data_description` for its explanation.

        show_hidden
          This parameter has no effect here. See :meth:`BaseDefinition.data_description` for its explanation.

        prefix
          Prefix for the indentation of the description.

        Returns
        -------
        additional_data_description

          An additional description of the values accepted by this configuration option, retrieved from the documentation type.
    """
    return self.type.additional_description(prefix)

  def validate(self, value, why='set'):
    if value is None:
       if self.required:
          raise ValueError(f"The value is required for {self.name}, cannot set it to None")
       return True
    if self.is_fixed and not np.array_equal(self.default_value, value):
       raise ValueError(f'The value of {self.name} is required to be {self.default_value}, cannot set it to {value}')
    self.type.validate(value, self.name, why=why)

  def convert_and_validate(self, value, why='set'):
    value = self.type.convert(value)
    self.validate(value, why)
    return value

  @property
  def value_name_format(self):
    return self.name_format

  @value_name_format.setter
  def value_name_format(self, value):
      self.name_format = value

  def __str__(self):
    out="<{}: {}>".format(self.name, str(self.type))
    val = self.get_value()
    if val is not None:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar(self, allow_dangerous=False):
    body = self.type.grammar(self.name)

    if self.is_fixed:
      def check_fixed(s, loc, x, body=body):
          if self.default_value.__class__ is np.ndarray:
             eq = np.array_equal(x[0], self.default_value)
          else:
             eq = x[0]==self.default_value
          if eq:
             return x
          message="The value of {} is {} and it should be {}".format(self.name, x[0], self.default_value)
          raise pp.ParseException(s,loc,message, body)
      body=body.copy().addParseAction(check_fixed)

    if allow_dangerous and hasattr(self, 'type_of_dangerous'):
        danger = pp.Forward()
        danger << self.type_of_dangerous.grammar(self.name + '_dangerous')
        danger.addParseAction(lambda x: DangerousValue(x[0], self.type_of_dangerous, False))
        body = body ^ danger

    if self.name_in_grammar:
       god = self.grammar_of_delimiter()
       body = god + body

    optional, df, _ = self.type.missing_value()
    if optional:
      body = pp.Optional(body).setParseAction( lambda x: x or df )
      nbody=''
    else:
      nbody=self.type.grammar_name()
      if self.name_in_grammar:
         nbody = str(god) + nbody

    out = self._tuple_with_my_name(body, has_value=self.type.has_value, is_numbered_array=self.is_numbered_array)
    out.setName(self.name + nbody)
    return out

  def get_value(self, option=None):
     """ Return the default or fixed value of this option.

     The function can accept the Option (which of the definition is): if the default value is given by callable,
     this argument is passed to it. (E.g. to set the default value using some properties obtained from the
     configuration objects.
     """
     if self.default_value is not None:
        if callable(self.default_value):
           return self.default_value(option)
        return self.default_value
     return None

  def write(self, file, value):
     """
     Write the option to the open file

     Parameters
     ----------
     file
      The file to write to.

     value
      The value to write. It can be instance of DangerousValue, in such case
      It's own type is used to write the value.
     """
     if self.type.has_value:
        missing, df, np = self.type.missing_value()

     def write(name, value):
         if isinstance(value, DangerousValue):
            type = value.value_type
            value = value()
         else:
            type = self.type

         if type.has_value:
             if value is None:
               value = self.get_value()
             if value is None:
               return False
             if np.__class__ is value.__class__ and np == value:
               return False
             write_value = not ( missing and df == value )
         else:
            value=None
            write_value=True

         file.write(self.prefix)
         if self.name_in_grammar:
            file.write(name)
            if write_value:
               file.write(self.name_value_delimiter)
               type.write(file, value)
         elif write_value:
            type.write(file, value)
         return True

     name = self.formated_name
     if self.is_numbered_array:
        out = False
        for i, val in value.items():
            out = write(name + (str(i) if i!='def' else ''), val) or out
        return out
     else:
        return write(name, value)

  def copy(self, **kwargs):
     default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
     default.update(kwargs)
     return self.__class__(**default)

  def remove(self, name):
     del self.section[name]
     return self

  def _generic_info(self):
      return f"Configuration value {self.name}"

  def can_be_repeated(self):
      return self.is_numbered_array

  def _get_copy_args(self)->Dict[str, str]:
       """
       Compute the dictionary that defines the attributes to create a copy of this object.

       Returns
       -------
       copy: Dict
          The returning dictionary has this structure:
          { name of the argument of the __init__ function : name of the object attribute }
       """
       out = super()._get_copy_args()
       if self.is_fixed:
           out['fixed_value'] = out['default_value']
       return out

  _copy_excluded_args = BaseDefinition._copy_excluded_args + ['fixed_value']

  def copy_value(self, value, all_values=False):
      """ Creates copy of the value

      Parameters
      ----------
      values
        The value to copy

      all_values
        Wheter, for a numbered array, a whole dict is supplied
      """
      if not all_values or not self.is_numbered_array:
          return self.type.copy_value(value)
      return { k:self.type.copy_value(v) for v in values }

def add_excluded_names_condition(element, names):
    """ Add the condition to the element, that
    its value is not any of given names """
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class ContainerDefinition(BaseDefinition):
    """ Base class for a definition (of contained data, format, etc)
    of either a whole configuration file
    (e.g. :class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` or
    e.g. :class:`Potential<ase2sprkkr.potentials.potentials.Potential>`) or
    its :class:`Section<ase2sprkkr.common.configuration_containers.Section>`.
    """

    force_order = False
    """ Force order of its members """

    value_name_format = None
    """ The (print) format, how the name is written """

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or OrderedDict()
        for value in args:
           items[value.name] = value
        return items

    def __init__(self, name, members=[], alternative_names=[], info=None, description=None,
                 is_optional=False, is_hidden=False, is_expert=False,
                 has_hidden_members=False, name_in_grammar=None, force_order=None,
                 write_alternative_name:bool=False, result_class=None):
       super().__init__(
           name = name,
           alternative_names = alternative_names,
           is_optional = is_optional,
           is_hidden = is_hidden,
           is_expert = is_expert,
           name_in_grammar = name_in_grammar,
           info = info,
           description = description,
           write_alternative_name = write_alternative_name,
           result_class = result_class
       )

       self.is_hidden = is_hidden
       if not isinstance(members, OrderedDict):
          members = self._dict_from_named_values(members)

       if self.value_name_format:
          for i in members.values():
              i.value_name_format = self.value_name_format
       self._members = members

       self.has_hidden_members = has_hidden_members
       if force_order is not None:
          self.force_order = force_order

    configuration_type_name = 'SECTION'
    """ Name of the container type in the runtime documentation """

    def data_description(self, verbose:Union[bool,str,int]=False, show_hidden:bool=False, prefix:str=''):
        """
        Return the runtime documentation for the configuration described by this object.

        Parameters
        ----------
        verbose
          If ``False``, only one line with the section name and basic info is returned.
          If ``True``, the items contained in the section are listed.
          If ``'all'``, add also detailed info about all descendants.
          If an ``int`` is given, print detailed informations about n levels. I.e. ``1`` is the same as ``True``

        show_hidden
          If False, do not show hidden members.

        prefix
          The string, with with each line will begin (commonly the spaces for the indentation).
        """
        def container_name():
            out = self.configuration_type_name
            if out:
               out+=' '
            out+=self.name
            return out

        out = f"{container_name()}"
        flags = []
        if self.force_order:
           flags.append('order fixed')
        if self.is_hidden:
           flags.append('hidden')
        if self.is_optional:
           flags.append('optional')
        if self.is_expert:
           flags.append('expert')
        if flags:
           flags = ', '.join(flags)
           out+=f" ({flags})"

        if verbose:
           if isinstance(verbose, int):
              verbose-=1
           else:
              verbose = verbose if verbose=='all' else False
           add = self.additional_data_description(verbose, show_hidden, prefix)
           if add:
              out+=f' contains:'
              under=prefix + "-"*len(out) + '\n'
              out=f"{prefix}{out}\n{under}{add}"
        return out

    def additional_data_description(self, verbose:Union[bool,str,int]=False, show_hidden=False, prefix:str=''):
        """
        Return the description (documentation for runtime) of the items in the container.

        Parameters
        ----------
        verbose
          If ``True``, include detailed description of the children.
          If ``'all'``, include even detailed description.
          If ``int`` is given, print detailed informations up to n levels.

        show_hidden
          If False, do not show hidden members.

        prefix
          The string, with with each line will begin (commonly the spaces for the indentation).
        """
        cprefix=prefix + self._description_indentation
        out = []

        def write(i):
           s = i.data_description(verbose, show_hidden, cprefix)
           if not '\n' in s:
              info = i.info(False)
              if info:
                 s = s + (' ' * (max(40 - len(s), 0) + 2)) + info
           else:
              s+='\n'
           out.append(s)

        expert = False
        for i in self:
            if i.is_hidden and not show_hidden:
               continue
            if not i.is_expert:
               write(i)
            else:
               expert=True

        if expert:
          out.append(f'\n{cprefix}Expert options:')
          out.append(cprefix +   '--------------')
          cprefix+=self._description_indentation

          for i in self:
              if i.is_expert:
                 if i.is_hidden and not show_hidden:
                    continue
                 write(i)

        return '\n'.join(out)

    def __iter__(self):
        return iter(self._members.values())

    def members(self):
        return self._members.values()

    def names(self):
        return self._members.keys()

    def __getitem__(self, key):
        return self._members[key]

    def __setitem__(self, key, value):
        self._members[key]=value

    def __contains__(self, key):
        return key in self._members

    def remove(self, name):
        del self._members[name]
        return self

    def copy(self, args=[], items=[], remove=[], defaults={}, **kwargs):
        """ Copy the section with the contained values modified by the arguments."""
        members = OrderedDict( ( (k,i.copy()) for k,i in self._members.items() ) )
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = members[i].type.convert(v)

        default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
        default.update(kwargs)
        default['members'] = members
        return self.__class__(**default)

    def create_object(self, container=None):
        return self.result_class(self, container)

    def all_member_names(self):
        return itertools.chain(
          (i.name for i in self),
          itertools.chain.from_iterable((
            i.alternative_names for i in self if i.alternative_names
          ))
        )

    def _values_grammar(self, allow_dangerous:bool=False, delimiter=None):
       if self.custom_class:
          custom_value = self.custom_member_grammar(self.all_member_names())
       else:
          custom_value = None
       delimiter = delimiter or self.grammar_of_delimiter()

       def grammars():
         """ This function iterates over the items of the container, joining all the without name_in_grammar with the previous ones. """

         def repeated_grammars():
             """ If the item can be repeated, do it here - we don't know, whether there is a fixed order in any way
             (e.g. the item is followed by the items without name in grammar)
             """
             for i in self._members.values():
                 if i.is_generated:
                     continue
                 g = i._grammar(allow_dangerous)
                 if i.can_be_repeated():
                     g = delimitedList(g, delimiter)
                 yield i,g

         it = iter(repeated_grammars())
         head_item, grammar_chain = next(it)

         for item, grammar in it:
             if item.name_in_grammar:
               yield head_item, grammar_chain
               head_item, grammar_chain = item, grammar
             else:
               add = delimiter + grammar
               if item.is_optional:
                  add = pp.Optional(add)
               grammar_chain = grammar_chain + add
         yield head_item, grammar_chain

       if self.force_order:

           init = pp.Empty()
           def set_loc(loc, toks):
               init.location = loc
           init.setParseAction(set_loc)

           first = pp.Empty().addCondition(lambda loc, toks: loc == init.location)
           if custom_value:
              cvs = pp.ZeroOrMore(custom_value + delimiter).setName('<custom...>')
              after = delimiter + cvs
           else:
              after = pp.Forward() << delimiter
           after.addCondition(lambda loc, toks: loc != init.location)
           inter = (first | after).setName('?')

           def sequence():
               for head,g in grammars():
                   g = inter + g
                   if head.is_optional:
                      g = pp.Optional(g)
                   yield g

           values  = pp.And([ i for i in sequence()])

           if custom_value:
              if not self._first_section_has_to_be_first():
                  values = cvs + values
              values += pp.ZeroOrMore(delimiter + custom_value)
           values = init + values

       else:
           it = grammars()
           #store the first fixed "chain of sections"
           first = self._first_section_has_to_be_first() and next(it)[1]
           #the rest has any order
           values = pp.MatchFirst([i for head,i in it])
           if custom_value:
               values |= custom_value
           values = delimitedList(values, delimiter)
           if first:
              values = first + pp.Optional(delimiter + values)

       values.setParseAction(lambda x: unique_dict(x.asList()))

       if self.validate:
          def _validate(s, loc, value):
              #just pass the dict to the validate function
              is_ok = self.validate(MergeDictAdaptor(value[0], self), 'parse')
              if is_ok is not True:
                if is_ok is None:
                   is_ok = f'Validation of parsed data of {self.name} section failed'
                raise pp.ParseException(s, loc, is_ok)
              return value
          values.addParseAction(_validate)

       return values

    def _grammar(self, allow_dangerous=False):
       delimiter = self.grammar_of_delimiter()
       values = self._values_grammar(allow_dangerous, delimiter)
       out = self._tuple_with_my_name(values, delimiter)
       out.setName(self.name)

       return out

    validate = None
    """ A function for validation of just the parsed result (not the user input) """

    @classmethod
    @cache
    def delimited_custom_value_grammar(cls):
        """ Return the grammar for the custom child with delimiter.
        The delimiter can delimite it either from the previous child or from the section name."""

        return cls.child_class.grammar_of_delimiter() + cls.custom_value_grammar()

    custom_name_characters = pp.alphanums + '_-()'
    """ Which characters can appears in an unknown child (value/section) name """

    @classmethod
    def custom_member_grammar(cls, value_names = []):
       """ Grammar for the custom - unknown - child """
       name = pp.Word(cls.custom_name_characters).setParseAction(lambda x: x[0].strip())
       add_excluded_names_condition(name, value_names)
       out = (name + cls.delimited_custom_value_grammar()).setParseAction(lambda x: tuple(x))
       out.setName(cls.custom_value_name)
       return out


    def _first_section_has_to_be_first(self):
       """ Has/ve the first child(s) in an unordered sequence fixed position? """
       return not self._members.first_item().name_in_grammar

    def parse_file(self, file, return_value_only=True, allow_dangerous=False):
       """ Parse the file, return the parsed data as dictionary """
       grammar = self.grammar(allow_dangerous)
       out = grammar.parseFile(file, **{ _parse_all_name: True } )
       return self.parse_return(out, return_value_only)

    def parse(self, string, whole_string=True, return_value_only=True, allow_dangerous=False):
       """ Parse the string, return the parsed data as dictionary """
       grammar = self.grammar(allow_dangerous)
       out = grammar.parseString(string, **{ _parse_all_name: whole_string } )
       return self.parse_return(out, return_value_only)

    def parse_return(self, val, return_value_only:bool=True):
        """ Clean up the parsed values (unpack then from unnecessary containers)

        Parameters
        ----------
        return_value_only
          Return only value, not name - value tuple

        """
        val = val[0]
        if return_value_only:
           val = val[1]
        return val

    async def parse_from_stream(self, stream, up_to, start=None, whole_string=True, return_value_only=True, allow_dangerous=False):
        """
        Parse string readed from asyncio stream.
        The stream is readed up to the given delimiter
        """

        result = await stream.readuntil(up_to)
        result = result[:-len(up_to)].decode('utf8')
        if start:
           result = start + result
        return self.parse(result, whole_string)

    def read_from_file(self, file, allow_dangerous=False, **kwargs):
        """ Read a configuration file and return the parsed Configuration object """
        out = self.result_class(definition = self, **kwargs)
        out.read_from_file(file, allow_dangerous=allow_dangerous)
        return out

    def validate(self, container, why:str='save'):
        return True

class SectionDefinition(ContainerDefinition):
   """ Base class for definition of the sections in Pot or InputParameters files.

       It just redefine a few properties/methods to values/behavior typical for the sections
   """

   result_class = Section

   @property
   def values(self):
       return self._members

   custom_value_name = 'CUSTOM_VALUE'
   """ Just the name that appears in the grammar, when it is printed."""

   @classmethod
   @cache
   def delimited_custom_value_grammar(cls):
        gt = cls.custom_class.grammar_type
        #here the child (Value) class delimiter should be used
        out = cls.child_class.grammar_of_delimiter() + gt.grammar()
        optional, df, _ = gt.missing_value()
        if optional:
           out = out | pp.Empty().setParseAction(lambda x: df)
        return out

   def _generic_info(self):
      return f"Configuration section {self.name}"


class ConfigurationRootDefinition(ContainerDefinition):
   """ From this class, the definition of the format of a whole configuration file should be derived.

   """

   name_in_grammar = False
   """ The configuration files has commonly no "name" of its content, they
   just contains their content.

   However, in some cases the name_in_grammar could be used for some kind of
   prefix in the file, however, it is better to use a fixed value for this purpose.
   """

   @classmethod
   def from_dict(cls, name, defs=None):
       """
       Create instance of the definition from a dictionary, creating
       the sections (and values) definitions recursively.
       """
       def gen(i):
           section = defs[i]
           if not isinstance(defs, SectionDefinition):
              section = cls.child_class(i, section)
           return section

       if defs is None:
          defs = name
          name = cls.__name__

       return cls(( gen(i) for i in defs))

   def __init__(self, name, members=[], **kwargs):
       if not members and not isinstance(name, str):
          members = name
          name = self.__class__.__name__
       super().__init__(name, members, **kwargs)

   @property
   def sections(self):
       return self._members

   custom_value_name = 'CUSTOM_SECTION'
   """ Just the name that appears in the grammar, when it is printed."""

   def _tuple_with_my_name(self, expr, delimiter=None):
       """ Do not create tuple (name, value) for the root class. """
       return expr

   def parse_return(self, val, return_value_only=True):
        """ Clean up the parsed values (unpack then from unnecessary containers)

            There is no name in the parsed results (see how
            ConfigurationRootDefinition._tuple_with_my_name is redefined).
        """
        val = val[0]
        return val

   def _grammar(self, allow_dangerous=False):
       """Returns the grammar to parse the configuration file.

       This method just tweaks the grammar (generated by the common container implementation) to ignore comments,
       so the comments would be ignored just once.
       """

       out=super()._grammar(allow_dangerous)
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out

class MergeDictAdaptor:
    """ This class returns a read-only dict-like class
    that merge values from a container and from the
    definition of a section """

    def __init__(self, values, definition):
        self.values = values
        self.definition = definition

    def __hasitem__(self, name):
        return self.values.__hasitem__(name) or \
               self.definition.__hasitem__(name)

    def __getitem__(self, name):
        try:
            return self.values[name]
        except KeyError:
            return self.definition[name].get_value()

    def __repr__(self):
        return f'Section {self.definition.name} with values {self.values}'
