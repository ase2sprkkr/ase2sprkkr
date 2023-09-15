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

from ..common.misc import dict_first_item
from ..common.grammar_types  import type_from_type, type_from_value, GrammarType, Array
from ..common.grammar  import delimitedList, end_of_file, generate_grammar
from .configuration_containers import Section
from .options import Option, Dummy, DangerousValue
from .decorators import cache, cached_class_property, cached_property

import numpy as np
import pyparsing as pp
import inspect
import itertools
import builtins
from io import StringIO
from typing import Dict, Union
import sys

#:This serves just for dealing with various pyparsing versions
_parse_all_name = 'parse_all' if \
  'parse_all' in inspect.getfullargspec(pp.Or.parseString).args \
  else 'parseAll'

def dict_from_parsed(values, allowed_duplicates):
    """ Create a dictionary from the arguments.
    From duplicate arguments create numpy arrays.
    Moreover, if there is key of type (a,b), it will be transformed to subdictionary.
    Such a keys do not allow duplicates.

    >>> dict_from_parsed( [ ('x', 'y'), (('a','b'), 1 ), (('a', 'c'), 2) ], [] )
    {'x': 'y', 'a': {'b': 1, 'c': 2}}
    >>> dict_from_parsed( [ ('x', 1), ('x', '2') ], [] ) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    pyparsing.exceptions.ParseException: There are non-unique keys: x
    >>> dict_from_parsed( [ ('x', 1), ('x', 2) ], ['x'] )
    {'x': array([1, 2])}
    """
    out = {}
    duplicates = []
    errors = []

    if isinstance(allowed_duplicates, list):
       _allowed_duplicates = lambda x: x in allowed_duplicates
    else:
       _allowed_duplicates = allowed_duplicates

    def add(out, k, key, check):
         if k in out:
            if key in duplicates:
               out[k].append(v)
            else:
               if not _allowed_duplicates(check):
                   return errors.append(check)
               out[k] = [ out[k], v ]
               duplicates.append(k)
         else:
            out[k] = v

    for k, v in values:
        if isinstance(k, tuple):
           if not k[0] in out:
              o = out[k[0]] = {}
           elif not isinstance(out[k[0]], dict):
              raise pp.ParseException(f"There is duplicate key {k[0]}, however, I should created "
                                "both arrays and dict from the data, which is not possible")
           else:
              o = out[k[0]]
           add(o, k[1], k, k[0])
        else:
           add(out, k, k, k)
    for k in duplicates:
        if isinstance(k, tuple):
            o=out[k[0]]
            k=k[1]
        else:
            o=out
        o[k] = np.array(out[k])
    if errors:
       errors = ','.join(errors)
       raise pp.ParseException(f"There are non-unique keys: {errors}")
    return out


class BaseDefinition:
  """ This class is a member of definition of configuration, that can be both
  real (holds a value or values) or just virtual.
  Virtual definitions do not store values, generates grammar e.g. using other members
  of the container etc...

  Parameters
  ----------
  name: str
      Name of the value/section

  is_optional: bool
      This element can be omited in the parsed file

  condition
     If defined, the condition
      - the condition.parse_condition() is invoked, when given grammar element
        should be parsed. If it is False, the element is skipped
      - the condition() is invoked, when the elements of the container is listed
        to hide the inactive members
  """

  def __init__(self, name, is_optional=False, condition=None):
       self.name = name
       self.is_optional=is_optional
       self.condition = condition
       self.grammar_hooks = []
       self.condition = condition

  def add_grammar_hook(self, hook):
       """ Added hooks process the grammar of the option/container.
       E.g. it is used when the number of readed lines should depend
       on the value of the option.
       """
       self.grammar_hooks.append(hook)

  def remove_grammar_hook(self, hook):
       """ Remove a grammar hook """
       self.grammar_hooks.remove(hook)

  def grammar(self, allow_dangerous:bool=False):
       """ Generate grammar with the correct settings of pyparsing global state.

       Parameters
       ----------
       allow_dangerous
        Allow dangerous values - i.e. values that do not fulfill the requirements
        for the given-option value (i.e. a type requirement or other constraints).
       """
       with generate_grammar():
         return self._grammar and self._grammar(allow_dangerous)

  @property
  def _grammar(self):
       """ Return the grammar. Descendants can redefine this method e.g. to allow
       not to generate the grammar at all.
       This method is implemented by property, that conditionaly returns the
       real method.

       Returns
          func: callable
          Function to generate grammar or None if no grammar is returned
       """
       return self._hooked_grammar

  def has_grammar(self):
       """ Returns, whether the definition generates a grammar or not. By default,
       it just check the self._grammar """
       return self._grammar

  def _hooked_grammar(self, allow_dangerous:bool=False, **kwargs):
       """ Generates grammar. Unlike :meth:`grammar`, it does not change the
       global pyparsing state to ensure that the generated grammar will handle
       whitespaces in a propper way. Unlike the :meth:`_create_grammar` method,
       which should contain implementation of the grammar creation, this function
       add the common functionality to the generated grammar (currently,
       just the grammar hooks)
       """
       out = self._create_grammar(allow_dangerous, **kwargs)
       return self._add_hooks_to_grammar(out)

  def _add_hooks_to_grammar(self, grammar):
       """ Add registered grammar hooks to a grammar """
       if self.grammar_hooks:
           for i in self.grammar_hooks:
               grammar=i(grammar)
       return grammar

  def added_to_container(self, container):
       """ Hook called, when the object is assigned to the container (currently from the container
       constructor) """
       self.container=container

  def info(self, *args, **kwargs):
       return "This object is not intended for a direct use."

  def description(self, *args, **kwargs):
       return "This object is not intended for a direct use."

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

  def copy(self, **kwargs):
     default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
     default.update(kwargs)
     return self.__class__(**default)

  def create_object(self, container=None):
     """ Creates Section/Option/.... object (whose properties I define) """
     return self.result_class(self, container)

  can_be_repeated = False
  """ If True, the item can be repeated in the parsed file.
  This behavior have currently the values with is_numbered_array property = True.
  This function is to be redefined in descendants
  """
  is_independent_on_the_predecessor = True

  _copy_excluded_args = ['container', 'grammar_hooks']

  def enrich(self, option):
        return option

  @property
  def output_definition(self):
        return self.__dict__.get('_output_definition', self)

  @output_definition.setter
  def output_definition(self, od):
        """ Set a special object for writing """
        self.__dict__['_output_definition'] = od


class RealItemDefinition(BaseDefinition):
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

   def __init__(self, name, written_name=None, alternative_names=None,
                is_optional=False, is_hidden=False, is_expert=False,
                name_in_grammar=None, info=None, description=None,
                write_alternative_name:bool=False,
                condition=None, write_condition=None,
                result_class=None
                ):
       """
       Parameters
       ----------
        name: str
          Name of the value/section

        written_name: str or None
          Name to write to the input file. Default None means use the name.

        alternative_names: str or [str]
          Alternative names that can denotes the value. If no written_name is given,
          the first alternative_names is used for the output. However, contrary to
          written_name, such way still allow to parse the name during parsing as
          the name of the value.

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

        condition
           If defined, the condition
            - the condition.parse_condition() is invoked, when given grammar element
              should be parsed. If it is False, the element is skipped
            - the condition() is invoked, when the elements of the container is listed
              to hide the inactive members

        result_class
           Redefine the class that holds data for this option/section.
       """
       super().__init__(name, is_optional, condition)
       self.written_name = written_name
       """ The name of the option/section """
       if isinstance(alternative_names, str):
          alternative_names = [ alternative_names ]
       self.alternative_names = alternative_names
       """ Alternative names of the option/section. The option/section can
       be "denoted" in the configuration file by either by its name or any
       of the alternative names.
       """
       self.is_expert = is_expert
       self.is_hidden = is_hidden
       """ Is it required part of configuration (or can it be ommited)? """
       self.write_alternative_name = write_alternative_name
       self.write_condition = write_condition or (lambda x: True)
       self.name_in_grammar = self.__class__.name_in_grammar \
                               if name_in_grammar is None else name_in_grammar
       self._info = info
       """ A short help text describing the content for the users. """
       self._description = description
       """ A longer help text describing the content for the users. """
       if result_class:
           self.result_class = result_class

   def all_names_in_grammar(self):
       if not self.name_in_grammar:
           return
       if self.written_name:
           yield self.written_name
       else:
           yield self.name
       if self.alternative_names:
           for i in self.alternative_names:
               yield i

   def allow_duplication(self):
       """ Can be the item repeated in the output file """
       return False

   @property
   def is_independent_on_the_predecessor(self):
       """ Some value have to be positioned just after their predecessor
       in the output.
       """
       return self.name_in_grammar

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

   def _grammar_of_name(self, is_numbered_array:bool=False):
        """
        Return grammar for the name (and possible alternative names etc.)

          Parameters
          ----------

          is_numbered_array
             If True, the resulting grammar is in the form
             NAME[index]
        """
        if self.name_in_grammar:
            names = self.all_names_in_grammar()
            keyword = pp.CaselessLiteral if is_numbered_array else pp.CaselessKeyword
            if self.do_not_skip_whitespaces_before_name:
               names = [ keyword(i).leaveWhitespace() for i in names ]
            else:
               names = [ keyword(i) for i in names ]
            if len(names) > 1:
                name = pp.Or(names)
                if self.do_not_skip_whitespaces_before_name:
                    names=names.leaveWhitespace()
            else:
                name = names[0]
            name.setParseAction(lambda x: self.name)
            if is_numbered_array:
               name += pp.Optional(pp.Word(pp.nums), default='def')
               name.setParseAction(lambda x: (x[0], 'def' if x[1]=='def' else int(x[1])) )
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)
        return name

   def _tuple_with_my_name(self, expr,
                           delimiter=None,
                           has_value:bool=True,
                           is_numbered_array:bool=False,
                           name_in_grammar=None):
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
        if self.name_in_grammar if name_in_grammar is None else name_in_grammar:
           name = self._grammar_of_name(is_numbered_array)
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

   _copy_excluded_args = BaseDefinition._copy_excluded_args + ['expert']

class ValueDefinition(RealItemDefinition):

  result_class = Option

  name_in_grammar = None
  is_generated = False

  def __init__(self, name, type=None, default_value=None,
               default_value_from_container=None,
               written_name=None, alternative_names=None,
               fixed_value=None, required=None, init_by_default=False,
               result_is_visible=False, info=None, description=None,
               is_stored=None, is_hidden=False, is_optional=None, is_expert=False,
               is_numbered_array:bool=False, is_repeated=False,
               is_always_added:bool=None,
               name_in_grammar=None, name_format=None, expert=None,
               write_alternative_name:bool=False,
               write_condition=None, condition=None,
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

    written_name: str
      Name of the configuration value in the input file

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

    init_by_default: bool
      If the value is not set, init it by default

    result_is_visible: bool
      If True, the result (see the result property of :class:`Option`) assigned to the
      option is visible if the value is get as its default value - in this mode the result
      can be used as a kind of default value, specific for given configuration object.

      If False, the result is visible to the user just using the result property
      and it should be some transformation of the value of the object (e.g.
      relative path of the given absolute path, id of assigned object etc.)

    is_stored: bool
      If True, the value is readed/writed from the output file.
      If None, set to False if the value is Generated.

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

    is_repeated
      If True, the real type of the Value is array of values of given type.
      The name-value pair is repeated for each value of the array.

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

    condition
       If defined, the condition
        - the condition.parse_condition() is invoked, when given grammar element
          should be parsed. If it is False, the element is skipped
        - the condition() is invoked, when the elements of the container is listed
          to hide the inactive members

    result_class
       Redefine the class that holds data for this option/section
    """
    if default_value_from_container:
       default_value = lambda o: default_value_from_container(o._container)
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

    self.init_by_default = init_by_default
    self.is_stored = not self.is_generated if is_stored is None else is_stored

    if fixed_value is None:
       self.is_fixed = False
    else:
       default_value = fixed_value
       self.is_fixed = True

    if default_value is None and not isinstance(type, (GrammarType, builtins.type)):
       self.type = type_from_value(type, type_map = self.type_from_type_map)
       self.default_value = None if isinstance(type,dict) else self.type.convert(type)
    else:
       self.type = type_from_type(type, type_map = self.type_from_type_map)
       if default_value is not None:
          if not callable(default_value):
             default_value = self.type.convert(default_value)
       self.default_value = default_value

    assert isinstance(self.type, GrammarType), "grammar_type (sprkkr.common.grammar_types.GrammarType descendat) required as a value type"
    self.type.used_in_definition(self)

    if is_repeated is True:
        self.is_repeated = self.type
        self.type = Array(self.type)
    else:
        self.is_repeated = is_repeated

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    if result_is_visible:
       self.default_value = lambda o: o._result if hasattr(o, '_result') else default_value

    if required is None:
       required = not is_expert and (not is_optional and default_value is None)
    self.required = required

    if is_optional is None:
       is_optional = required is False

    super().__init__(
         name = name,
         written_name = written_name,
         alternative_names = alternative_names,
         is_optional = is_optional,
         is_hidden = is_hidden,
         is_expert = is_expert,
         name_in_grammar = name_in_grammar,
         info=info,
         description = description,
         write_alternative_name = write_alternative_name,
         write_condition = write_condition,
         condition = condition,
         result_class = result_class,
    )

    if self.name_in_grammar is None:
        self.name_in_grammar = self.type.name_in_grammar

    self.is_numbered_array = is_numbered_array
    if is_numbered_array and not self.name_in_grammar:
       raise ValueException('Numbered_array value type has to have its name in the grammar')

    self.name_format = name_format

  configuration_type_name = 'OPTION'

  type_from_type_map = {}
  """ Redefine this in descendants, if you need to create different types that the defaults to be
  'guessed' from the default values """

  def allow_duplication(self):
       """ Can be the item repeated in the output file """
       return self.is_repeated

  @property
  def is_independent_on_the_predecessor(self):
      """ Some value have to be positioned just after their predecessor
      in the output.
      """
      return self.name_in_grammar or self.type.is_independent_on_the_predecessor

  def enrich(self, option):
      """ The Option can be enriched by the definition, e.g. the docsting can be extended. """
      self.type.enrich(option)

  @property
  def formated_name(self):
    if self.written_name:
       name = self.written_name
    else:
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
          This parameter has no effect here. See :meth:`RealItemDefinition.data_description` for its explanation.

        show_hidden
          This parameter has no effect here. See :meth:`RealItemDefinition.data_description` for its explanation.

        prefix
          Prefix for the indentation of the description.

        Returns
        -------
        additional_data_description

          An additional description of the values accepted by this configuration option, retrieved from the documentation type.
    """
    return self.type.additional_description(prefix)

  def added_to_container(self, container):
       """ Hook called, when the object is assigned to the container (currently from the container
       constructor) """
       self.container = container
       self.type.added_to_container(container)

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
    try:
      val = self.get_value()
    except Exception:
      val = "<ERRORNEOUS VALUE>"
    if val is not None:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar_of_value(self, delimiter, allow_dangerous=False):
    """ Return grammar for the (possible optional) value pair """
    type = self.is_repeated or self.type
    body = type.grammar(self.name)

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

    if delimiter:
       body = delimiter + body

    optional, df, _ = type.missing_value()
    if optional:
      body = pp.Optional(body).setParseAction( lambda x: x or df )
    return body

  @property
  def _grammar(self):
     if not self.is_stored:
        return None
     return self._hooked_grammar

  def _create_grammar(self, allow_dangerous=False, name_in_grammar=None, name_value_delimiter=None, original=False):
    """ Return a grammar for the name-value pair """
    if self.output_definition is not self and not original:
        g = self.output_definition._grammar
        return g and g(allow_dangerous)

    name_in_grammar = name_in_grammar if name_in_grammar is not None else self.name_in_grammar
    if (name_value_delimiter is None and name_in_grammar) or name_value_delimiter is True:
       name_value_delimiter=self.grammar_of_delimiter

    body = self._grammar_of_value(name_value_delimiter, allow_dangerous)


    if name_in_grammar:
       nbody = self.formated_name.strip()
    else:
       nbody = ''

    if not self.type.missing_value()[0]:
        if nbody:
            nbody +=str(name_value_delimiter) or ' '
        nbody+=self.type.grammar_name()

    out = self._tuple_with_my_name(body, has_value=self.type.has_value,
                                         is_numbered_array=self.is_numbered_array,
                                         name_in_grammar=name_in_grammar)
    out.setName(nbody)
    if self.is_repeated:
        out = out + pp.ZeroOrMore(self.container.grammar_of_delimiter + out)
    return out

  def get_value(self, option=None):
     """ Return the default or fixed value of this option.

     The function can accept the Option (which of the definition is): if the default value is given by callable,
     this argument is passed to it. (E.g. to set the default value using some properties obtained from the
     configuration objects.
     """
     if self.default_value is not None:
        if callable(self.default_value):
           return self.type.convert(self.default_value(option))
        return self.default_value
     return None

  def _save_to_file(self, file, option, always=False):
      value, write = option._written_value(always)
      if write:
          return self.write(file, value)
      else:
          return

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
     def write(name, value):
         if self.name_in_grammar:
            self.write_name(file, name)
            return self.write_value(file, value, self.name_value_delimiter)
         else:
            return self.write_value(file, value, self.prefix)

     name = self.formated_name
     if self.is_numbered_array or self.is_repeated:

        if self.is_numbered_array:
            written = ( (name + (str(i) if i!='def' else ''), v) for i, v in value.items() )
        else:
            written = ( (name, v) for v in value )
        out = False
        for mname, val in written:
            if out:
               file.write(self.container.delimiter)
            out = write(mname, val) or out
        return out

        out = False
     else:
        return write(name, value)

  def write_value(self, file, value, delimiter=''):
     """ Write given value and a given delimiter before the value to a file """
     if isinstance(value, DangerousValue):
        type = value.value_type
        value = value()
     else:
        type = self.is_repeated or self.type

     if type.has_value:
         if value is None:
           value = self.get_value()
         if value is None:
           return False
         missing, df, _ = self.type.missing_value()
         write_value = not ( missing and df == value )
     else:
        value=None
        write_value=True

     if write_value:
         file.write(delimiter)
         type.write(file, value)
         return True
     else:
         return False

  def write_name(self, file, name):
      file.write(self.prefix)
      file.write(name)

  def remove(self, name):
     del self.section[name]
     return self

  def _generic_info(self):
      return f"Configuration value {self.name}"

  @property
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

  _copy_excluded_args = RealItemDefinition._copy_excluded_args + ['fixed_value', 'result_is_visible', 'default_value_from_container']

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

  @property
  def is_independent_on_the_predecessor(self):
      return self.name_in_grammar or self.type.is_independent_on_the_predecessor

def add_excluded_names_condition(element, names):
    """ Add the condition to the element, that
    its value is not any of given names """
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class ContainerDefinition(RealItemDefinition):
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

    write_last_delimiter = True

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or {}
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
       if not isinstance(members, dict):
          members = self._dict_from_named_values(members)

       if self.value_name_format:
          for i in members.values():
              i.value_name_format = self.value_name_format
       self._members = members
       for i in self._members.values():
           i.added_to_container(self)

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
        members = dict( ( (k,i.copy()) for k,i in self._members.items() ) )
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = members[i].type.convert(v)

        default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
        default.update(kwargs)
        default['members'] = members
        return self.__class__(**default)

    def all_member_names(self):
        return itertools.chain.from_iterable(
            i.all_names_in_grammar() for i in self
        )

    def _grammar_of_values(self, allow_dangerous:bool=False, delimiter=None):
       if self.custom_class:
          custom_value = self.custom_member_grammar(self.all_member_names())
       else:
          custom_value = None
       delimiter = delimiter or self.grammar_of_delimiter

       def repeated_grammars():
           """ If the item can be repeated, do it here - we don't know, whether there is a fixed order in any way
           (e.g. the item is followed by the items without name in grammar)
           """
           for i in self._members.values():
               g = i._grammar and i._grammar(allow_dangerous)
               if not g:
                   continue
               if i.can_be_repeated:
                   g = delimitedList(g, delimiter)
               yield i,g

       def grammars():
         """ This function iterates over the items of the container, joining all the without name_in_grammar with the previous ones. """
         it = iter(repeated_grammars())
         head_item, grammar_chain = next(it)

         for item, grammar in it:
             if item.is_independent_on_the_predecessor:
               yield head_item, grammar_chain
               head_item, grammar_chain = item, grammar
             else:
               add = delimiter + grammar
               if item.is_optional:
                  add = pp.Optional(add)
               if item.condition and self.force_order:
                  add = item.condition.prepare_grammar(item, add)
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
           inter_cvs = (first | after).setName(f'<?DELIM>')
           inter = (first | delimiter.copy().addCondition(lambda loc, toks: loc != init.location))

           def sequence():
               for head,g in repeated_grammars():
                   if head.is_independent_on_the_predecessor:
                      delim = inter_cvs
                   else:
                      delim = inter
                   g = delim + g
                   if head.is_optional:
                      g = pp.Optional(g)
                   if head.condition:
                      yield head.condition.prepare_grammar(head, g)
                   else:
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

       values.setParseAction(lambda x: dict_from_parsed(x.asList(), self._allow_duplicates_of))

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

    def _allow_duplicates_of(self, name):
        """ Can a given element (identified by name) have more values in the parsed results?
            (However, not all definitions have to specify allow_duplicates, just the ones
            that have a value). For the others, this function raises an error.
        """
        return self[name].allow_duplication()

    def _create_grammar(self, allow_dangerous=False):
       delimiter = self.grammar_of_delimiter
       values = self._grammar_of_values(allow_dangerous, delimiter)
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

        return cls.child_class.grammar_of_delimiter + cls.custom_value_grammar()

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
       return not dict_first_item(self._members).is_independent_on_the_predecessor

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

    def read_from_string(self, string, allow_dangerous=False, **kwargs):
        return self.read_from_file(StringIO(string), allow_dangerous, **kwargs)

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
        out = cls.child_class.grammar_of_delimiter + gt.grammar()
        optional, df, _ = gt.missing_value()
        if optional:
           out = out | pp.Empty().setParseAction(lambda x: df)
        return out

   def _generic_info(self):
      return f"Configuration section {self.name}"


class ConfigurationRootDefinition(ContainerDefinition):
   """ From this class, the definition of the format of a whole configuration file should be derived.

   """
   write_last_delimiter = False
   """ Do not print additional newline after the last section """

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

   def _create_grammar(self, allow_dangerous=False):
       """Returns the grammar to parse the configuration file.

       This method just tweaks the grammar (generated by the common container implementation) to ignore comments,
       so the comments would be ignored just once.
       """
       out=super()._create_grammar(allow_dangerous)
       out=self.add_ignored(out)
       return out

   def add_ignored(self, grammar):
       grammar = grammar.ignore("#" + pp.restOfLine + pp.LineEnd())
       return grammar

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

class VirtualDefinition(BaseDefinition):
     """ Base class for a definition, that do not have value, just control
     the flow of the parsing """
     is_hidden = True
     counter = 1

     def create_object(self, container=None):
         return Dummy(self, container)

     def __init__(self, name=None, template=None, condition=None):
         if not name:
             if not template:
                 template=self.__class__.__name__.upper()
             name = f"_{template}_{VirtualDefinition.counter}"
             VirtualDefinition.counter+=1
         super().__init__(name, condition)

     def all_names_in_grammar(self):
         return iter(())

class ControlDefinition(VirtualDefinition):
     """ Control definitions has no grammar, they just modify the other items of the container """

     _grammar = None
     """ This object does not generate a grammar """

class Stub(VirtualDefinition):
    """ Item that allows to reuse existing item on the other place
    e.g. in another branch of Switch """

    def __init__(self, item, name=None, condition=None):
        super().__init__(name=None, template=f'STUB FOR {name}', condition=None)
        self.item=item
        self.condition=None

    def _save_to_file(self, file, value, always=True):
         item = value._container[self.item]
         if self.condition and not self.condition(value):
             return
         return item._save_to_file(file, always=True)

    def _create_grammar(self, allow_dangerous=False, **kwargs):
         item = self.container[self.item]
         out=item._grammar(allow_dangerous)
         return pp.Forward() << out

class Ignored:
    """ Output definition for an ignored option.
        Output definition can override the standard definition and
        set a special way how the item is read/writen:
        such option is not readed/writed at all... or is readed/writed
        by an another option, see :meth:`gather`
    """

    @cached_class_property
    def singleton(cls):
        return cls()

    _grammar = None

    def has_grammar(self):
        return False

    def _save_to_file(self, file, value, always=False):
        return


class Gather:
    """
    Output definition for the element of grammar, that reads besides himself
    other grammar elements, such that their names goes first and then the
    values go. See gather."""

    def __init__(self, *items, name_delimiter=' ', value_delimiter='\t', value_delimiter_grammar=''):
        if len(items) == 0:
           raise ValueError('Gather requires at least one value')
        self.items = items
        self.name_delimiter = name_delimiter
        self.value_delimiter = value_delimiter
        self.value_delimiter_grammar = value_delimiter_grammar

    def _grammar(self, allow_dangerous=False, **kwargs):
        names = [ i._grammar_of_name() for i in self.items if i.name_in_grammar ]
        ln = len(names)
        delimiter = bool(names)

        def values():
            nonlocal delimiter
            for i in self.items:
                yield i._grammar(allow_dangerous,
                                 name_value_delimiter=delimiter,
                                 name_in_grammar=False,
                                 original=True)
                delimiter=self.value_delimiter_grammar
        names.extend(values())
        out=pp.And(names)

        def discard_names(x):
            x = x.asList()
            return x[ln:]

        out.setParseAction(discard_names)
        return out

    def _save_to_file(self, file, value, always=False):
        names = self.name_delimiter.join(i.formated_name for i in self.items if i.name_in_grammar)
        if names:
            value._definition.write_name(file, names)
            delimiter = self.items[0].name_value_delimiter
        else:
            delimiter = ''

        def write(i, delimiter):
            val,write = i._written_value()
            if write:
                if not i._definition.write_value(file, val, delimiter):
                    raise NotImplemented('Gathered values have to be always written')
            else:
                raise NotImplemented('Gathered value names have to be always written')
            delimiter = self.value_delimiter

        write(value, delimiter)
        for i in self.items[1:]:
          write(value._container[i.name], self.value_delimiter)
        return True

    def copy_for_value(self, new_value):
        out = copy.copy()
        parent = new_value._container
        self.values = [ parent[i.name] for i in self.values ]

def gather(first, *members):
    """ Modify the given option definitions, that they appears
    in the output file in the form::

        <NAME> <NAME 2> <NAME 3> ... = <VALUE 1> <VALUE 2> <VALUE 3> ...
    """
    first.output_definition = Gather(first, *members)
    for i in members:
        i.output_definition = Ignored.singleton

    return (first, ) + members


def switch(item, values, name=None):
    switch = Switch(item, values, name)
    return (switch, ) + tuple(switch.all_values())

class Switch(ControlDefinition):
   """ Items of this class can control, which elements of grammar will be active and which not """

   create_object = None

   with generate_grammar():
       empty = pp.Empty()

   def __init__(self, item, values, name=None, template=None):
       """
       Parameters
       ----------
       item
          The name of Option, whose value determine the active elements
       values
          Dictionary, with the possible values of the item in the keys and
          the active elements in the values.

          example::

            V('TYPE', Keyword('SCALAR', 'COMPLEX'),
            Switch(
                {'SCALAR' : V('SCALAR' : int),
               'COMPLEX': V('COMPLEX', complex) }

          desribes both the following files::

            TYPE=SCALAR
            SCALAR=1

          and::
            TYPE=COMPLEX
            COMPLEX=1e5 7e5

       name
          Not needed to be supplied, it can be autogenerated.
       """

       self.item = item
       used = set()

       def create(n):
           nonlocal used
           if n.name in used:
               return Stub(n.name, condition=self)
           else:
               used.add(n.name)
               return n

       def convert(v):
           nonlocal used
           if isinstance(v, dict):
                return { i: create(v) for i,v in v.items() }
           else:
                if not isinstance(v, (tuple, list)):
                   v=[v]
                return { i.name: create(i) for i in v }
       self.values = { k : convert(v) if isinstance(v, (list, tuple)) else { v.name: v} for k,v in values.items() }
       self.container = None
       self.copied = False
       if template is None:
          template='_SWITCH_{item}'
       super().__init__(name, template)

   def copy(self):
       raise NotImplemented

   def __call__(self, option):
       return option._definition in self.values[option._container[self.item]()].values()

   def all_values(self):
       return itertools.chain.from_iterable( (i.values() for i in self.values.values()) )

   def item_hook(self, grammar):
       if not self.container.force_order:
          raise NotImplemented("Switch for custom-order containers have not yet been implemented")
       def item_value(value):
           ok = self.values.get(value, {})
           for i in ok.values():
               tpl = grammar._prepared.get(i.name, None)
               if tpl:
                  tpl[1] << tpl[0]
                  tpl[1].setName(f"<IF True THEN {str(tpl[0])}>")
               elif i.output_definition.has_grammar():
                  raise KeyError(f"In Switch, the item {i.name} for case {value} was not prepared")

           #no.setParseAction(lambda x: breakpoint() or x)
           for i in set(self.all_values()).difference(ok.values()):
               tpl = grammar._prepared.get(i.name, None)
               if tpl:
                  tpl[1].setName(f"<IF False THEN {str(tpl[0])}>")
                  tpl[1] << self.empty
               elif i.output_definition.has_grammar():
                  raise KeyError(f"In Switch, the item {i.name} for case {value} was not prepared")
       grammar._prepared = {}
       #if getattr(sys, "hhh", None): breakpoint()
       self.grammar = grammar
       return grammar.addParseAction(lambda x: item_value(x[0][1]) and x)

   def prepare_grammar(self, definition, grammar):
       f = pp.Forward()
       #if getattr(sys, "hhh", None): breakpoint()
       #old pyparsing compatibility
       if hasattr(f, 'set_name'):
           f.set_name(f"<IF {self.item} THEN {grammar.name}>")
       self.grammar._prepared[definition.name] = (grammar, f)
       return f

   def remove_from_container(self):
       if self.container:
           self.container[self.item].remove_grammar_hook(self.item_hook)
           for i in self.values.values():
               for j in i.values():
                   j.condition = None

   def added_to_container(self, container):
       self.remove_from_container()
       if container:
           if self.copied:
              self.copied = False
              self.values = { k:{ n : container[n] for n in v } for k,v in self.values }

           container[self.item].add_grammar_hook(self.item_hook)
           for i in self.values.values():
                for j in i.values():
                    j.condition = self
       super().added_to_container(container)

   def __del__(self):
       self.remove_from_container()


class SeparatorDefinition(VirtualDefinition):
    """ Basic class for separators """
    def __repr__(self):
        return "<SEPARATOR>"

    def __init__(self):
        super().__init__(template='SEPARATOR')

    def _create_grammar(self, allow_dangerous=False):
        return pp.Suppress(self.separator_type.grammar())

    def _save_to_file(self, file, value, always=False):
        if not always:
            if self.condition and self.condition(container):
                    return False
            self.separator_type.write(file, None)
            return True
import sys
