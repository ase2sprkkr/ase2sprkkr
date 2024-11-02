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

import pyparsing as pp
import inspect
from typing import Dict, Union
import itertools
from . import backward_compatibility  # NOQA
from enum import Enum, nonmember
import numpy as np

from .warnings import DataValidityWarning
from .options import Dummy, DummyStub
from .decorators import cached_class_property
from .grammar import generate_grammar
from .grammar_types.basic import Separator
from .parsing_results import Key, ArrayKey, DictKey, RepeatedKey, DefDictKey, IgnoredKey


class NotAllowed:
    """ Object of this class can be returned as a result of a :func:ConfigurationDefinition.condition """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

    def __bool__(self):
        """ When returned, the value is NOT ALLOWED"""
        return False


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
  validate_section = None
  """ Can be redefined for validating whole section """

  class Repeated(Enum):

      @nonmember
      class Type(Enum):
          """ Type of repetition. """
          NO = None
          # No repetition allowed
          ARRAY = np.ndarray
          # The item can be repeated, the result is an (dense) array of values
          LIST = list
          # The item can be repeated, the result is an (dense) array of values
          DICT = dict
          # The item can be repeated, the result is a (sparse) dict of values

          def __bool__(self):
              return self.value

      @nonmember
      class Numbering(Enum):
          """ Is the output numbered or not? """
          NO = 0
          # No numbering, just repeating (if allowed)
          YES = 1
          # Yes, number
          WITH_DEFAULT = 2
          # Numbered output, default value without number

          def __bool__(self):
              return self.value > 0

          @property
          def has_default(self):
              return self.value == 2

      def __init__(self, type:Type, key_type:Union[Key, callable]=Key.NONE,
                   is_numbered:Numbering=Numbering.NO,
                   has_header:bool=True):
          """
          Params
          ------
          type
            Type of the repetition (and the resulting storage)

          key_type
            If grammar emits a special type for keys (name of the option), pass it here.
            Such special keys takes care of aggregating the values into the propper storage (dict, list,...)

          is_numbered
            Whether the repeated occurences are numbered or not and how. Numbered
            means, that there can be values in the form
            NAME1=... NAME2=.... etc...
          """

          self.type = type
          self.key_type = key_type
          self.is_numbered = is_numbered
          self.has_header = True

      def __bool__(self):
          return self.type != self.Type.NO

      @classmethod
      def create(cls, val, grammar_type=None):
          if val is False:
              return cls.NO
          if val is True:
              if hasattr(grammar_type, 'is_numpy_array') and grammar_type.is_numpy_array:
                  return cls.ARRAY
              else:
                  return cls.REPEATED
          if isinstance(val, str):
              return cls[val]
          return val

      @property
      def is_array(self):
          return self.type in (self.Type.ARRAY, self.Type.LIST)

      @property
      def is_dict(self):
          return self.type == self.Type.DICT

      NO = (Type.NO)
      # No repetition at all
      IGNORED = (Type.NO, IgnoredKey)
      # The value is completely ignored
      REPEATED = (Type.LIST, RepeatedKey)
      # Values can be repeated, the result is array of values
      ARRAY = (Type.ARRAY, RepeatedKey)
      # Values can be repeated, the result is array of values
      REPEATED_SECTION = (Type.ARRAY, Key.NONE, Numbering.NO, False)
      # Repeated sections has no header, and their repetition solves themselves
      NUMBERED = (Type.ARRAY, ArrayKey, Numbering.YES)
      # Values are given in form "{NAME}{INDEX}", the result is array of values
      DICT = (Type.DICT, DictKey, Numbering.YES)
      # Values are given in form "{NAME}{INDEX}", the result is dict of values
      DEFAULTDICT = (Type.DICT, DefDictKey, Numbering.WITH_DEFAULT)
      # Same as NUMBERED_DICT, with non-numbered item possible, serving as default value

#
#
  def __init__(self, name, is_optional=False, condition=None):
       self.name = name
       self.name_lcase = name.lower()
       self.is_optional=is_optional
       self.condition = condition
       self.grammar_hooks = []
       self.condition = condition
       self.container = None

  is_repeated = Repeated.NO
  """ By default, the configuration items are not repeated """

  @property
  def real_name(self):
      return self.name

  def has_name(self, name, lower_case=False):
       return name == ( self.name_lcase if lower_case else self.name )

  def allowed(self, container):
        return not self.condition or self.condition(self, container)

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
       """ Add registered grammar hooks to a grammar.

       Developer note: the hooks are added even by grammar types, so it
       can not accept definition as argument, if you need it, use the closure
       of the hook to "pass" the argument.
       """
       if self.grammar_hooks:
           for i in self.grammar_hooks:
               grammar=i(grammar)
       return grammar

  def added_to_container(self, container):
       """ Hook called, when the object is assigned to the container (currently from the container
       constructor) """
       self.container=container

  def get_path(self):
      out = self.name
      if self.container and self.container.container:
          return self.container.get_path() + "." + out
      return out

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
          self.__class__._copy_args = { v: '_' + v if '_' + v in self.__dict__ else v
                                       for v in args if v not in self._copy_excluded_args }
       return self.__class__._copy_args

  def copy(self, **kwargs):
     default = { k: getattr(self, v) for k,v in self._get_copy_args().items() }
     default.update(kwargs)
     return self.__class__(**default)

  def create_object(self, container=None):
     """ Creates Section/Option/.... object (whose properties I define) """
     return self.result_class(self, container)

  can_be_repeated = False
  """ If True, the item can be repeated in the parsed file. The results will
  appear multiple times in the resulting dictionary after the parse.
  This attribute/property function is to be redefined in descendants.
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

  def accept_value(self, value) -> bool:
      """
      Return, whether a given value is accepted by this type of object.
      Accepting does not mean that the value can be validated. It is a
      base check, that a type of the value is suitable to be used here.
      If it is not -- e.g. single value is assigned to a section -- then
      the container can try to find another object (e.g. contained value
      of the same name) where the value is accepted.
      """
      return True


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
   info_in_data_description = False

   def __init__(self, name, written_name=None, alternative_names=None,
                is_optional=False, is_hidden=False, is_expert=False,
                name_in_grammar=None, name_format=None,
                info=None, description=None,
                write_alternative_name:bool=False,
                name_regex=None,
                condition=None, write_condition=None,
                result_class=None, warning_condition=None,
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

        name_format: str or None
          The way how the name is written

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
            - the condition.parse_condition() is invoked, when a given grammar element
              should be parsed. If it is False, the element is skipped
            - the condition() is invoked, when the elements of the container is listed
              to hide the inactive members

        result_class
           Redefine the class that holds data for this option/section.

        warning_condition
           If this lambda returns a non-none during validation, a warning will be issued.
       """
       super().__init__(name, is_optional, condition)
       self.written_name = written_name
       """ The name of the option/section """
       if isinstance(alternative_names, str):
           alternative_names = [ alternative_names ]
       self.alternative_names = alternative_names
       if alternative_names:
           self.alternative_names_lcase = [ i.lower() for i in alternative_names ]
       else:
           self.alternative_names_lcase = self.alternative_names
       self.name_regex = name_regex
       """ Alternative names of the option/section. The option/section can
       be "denoted" in the configuration file by either by its name or any
       of the alternative names.
       """
       self.is_expert = is_expert
       self.is_hidden = is_hidden
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
       self.warning_condition = None
       self.name_format = name_format

   @property
   def formated_name(self):
      if self.written_name:
         name = self.written_name
      else:
         name = next(iter(self.alternative_names)) if self.write_alternative_name else self.name
      if self.name_format:
         return "{:{}}".format(name, self.name_format)
      return name

   def has_name(self, name, lower_case=False):
       if super().has_name(name, lower_case):
           return True
       if self.alternative_names:
         if name in (self.alternative_names_lcase if lower_case else self.alternative_names):
           return True
       return False

   def validate_warning(self, value):
       if self.warning_condition:
          out = self.warning_condition(value)
          if out is not None:
              DataValidityWarning.warn(out)

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
       if not out:
          if generic:
             out = self._generic_info()
          else:
             out = ''
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

       # out.append(f"{prefix}Data description\n"
       #            f"{prefix}----------------")

       out.append(self.data_description(verbose, show_hidden, prefix))
       if self._description:
          out.append('\n')
          out.append(prefix + self._description.replace('\n', '\n' + prefix))
       return '\n'.join(out)

   def _grammar_of_name(self):
        """
        Return grammar for the name (and possible alternative names etc.)
        """
        if self.name_in_grammar:
            if self.name_regex:
                reg = pp.Regex(self.name_regex)
                if self.do_not_skip_whitespaces_before_name:
                   reg.leaveWhitespace()
                names=[pp.Regex(self.name_regex)]
            else:
                names = self.all_names_in_grammar()
                keyword = pp.CaselessLiteral if self.is_repeated.is_numbered else pp.CaselessKeyword
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
            if self.is_repeated:
                if self.is_repeated.is_numbered:
                    idx = pp.Word(pp.nums)
                    if self.is_repeated.is_numbered.has_default:
                      name += pp.Optional(idx, default='def')
                    else:
                      name += idx
                    name+= pp.WordEnd(pp.alphanums + "_")
            name.setParseAction(lambda x: self.is_repeated.key_type(self.name, *x.asList()[1:]))
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)

        return name

   def _tuple_with_my_name(self, expr,
                           delimiter=None,
                           has_value:bool=True,
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
        """
        if self.name_in_grammar if name_in_grammar is None else name_in_grammar:
           name = self._grammar_of_name()
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


class VirtualDefinition(BaseDefinition):
     """ Base class for a definition, that do not have value, just control
     the flow of the parsing """
     is_hidden = True
     counter = 1
     item_type = 'control item'

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
        super().__init__(name=None, template=f'STUB_FOR_{name or item}', condition=None)
        self.item=item
        self.condition=None

    def create_object(self, container=None):
         return DummyStub(self, container)

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter=''):
         item = value._container[self.item]
         if not always and not self.allowed(value._container):
             return False
         return item._save_to_file(file, always=True, name_in_grammar=name_in_grammar, delimiter=delimiter)

    def _create_grammar(self, allow_dangerous=False, **kwargs):
         item = self.container[self.item]
         out=item._grammar(allow_dangerous)
         return pp.Forward() << out

    @property
    def real_name(self):
        return self.item


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

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter=''):
        return False


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

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter=''):
        if name_in_grammar is not False:
            names = self.name_delimiter.join(i.formated_name for i in self.items if i.name_in_grammar)
        else:
            names = None
        if names:
            value._definition.write_name(file, names, delimiter)
            delimiter = self.items[0].name_value_delimiter
        else:
            delimiter = ''

        def write(i):
            val,write = i._written_value()
            if write:
                if not i._definition.write_value(file, val, delimiter):
                    raise NotImplementedError('Gathered values have to be always written')
            else:
                raise NotImplementedError('Gathered values have to be always written')

        write(value)
        delimiter = self.value_delimiter
        for i in self.items[1:]:
             write(value._container[i.name])
        return True


def gather(first, *members):
    """ Modify the given option definitions, that they appears
    in the output file in the form::

        <NAME> <NAME 2> <NAME 3> ... = <VALUE 1> <VALUE 2> <VALUE 3> ...
    """
    first.output_definition = Gather(first, *members)
    for i in members:
        i.output_definition = Ignored.singleton

    return (first, ) + members


def switch(item, values, condition=lambda x:x, name=None):
    switch = Switch(item, values, condition, name)
    return (switch, ) + tuple(switch.all_values())


def if_defined(item, values, not_values=[], name=None):
    if name is None:
        name = "IF_{item}"
    return switch(item, { None: not_values, True: values}, lambda x: None if x is None else True, name)


def if_not_defined(item, values, not_values=[], name=None):
    if name is None:
        name = "IF_{item}"
    return if_defined(item, not_values, values, name)


class Switch(ControlDefinition):
   """ Items of this class can control, which elements of grammar will be active and which not """

   create_object = None

   with generate_grammar():
       empty = pp.Empty()
       empty.__xx = 1

   def __init__(self, item, values, condition=lambda x:x, name=None, template=None):
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

       used = set()
       self.item = item
       self.fn = condition

       def create(n):
           nonlocal used
           if n.name in used:
               return Stub(n.name, condition=self)
           else:
               used.add(n.name)
               return n

       def convert(v):
           if isinstance(v, dict):
                return { i: create(v) for i,v in v.items() }
           else:
                if not isinstance(v, (tuple, list)):
                   v=[v]
                return { i.name: create(i) for i in v }
       self.values = { k : convert(v) if isinstance(v, (list, tuple)) else { v.name: v} for k,v in values.items() }
       self.container = None
       self.copied = False

       if template is None and name is None:
           name=f'_SWITCH_{item}'
       super().__init__(name, template)

   def copy(self):
       raise NotImplementedError

   def __call__(self, definition, data):
       return definition in self.items_for(self.switch_value(data))

   def items_for(self, value):
       return self.values[value].values()

   def all_values(self):
       return itertools.chain.from_iterable( (i.values() for i in self.values.values()) )

   def prepare_grammar(self, definition, grammar):
       """ Childs will have either the original grammar, or empty - set by _choose method.
       So replace the grammar by the forward.
       """
       if not self.container.force_order:
            return grammar
       f = pp.Forward()
       if definition.name in self.values.get(None, {}):
           f <<= grammar
           f.setName(f"<IF* ndef THEN {str(grammar)}>")
       else:
           f <<= self.empty
           f.setName(f"<IF* False THEN {str(grammar)}")
       self._grammars[definition.name] = (grammar, f)
       return f

   def item_hook(self, grammar):
       grammars = self._grammars = {}

       def choose(value):
           ok = self.values.get(value, {})
           for i in ok.values():
               tpl = grammars.get(i.name, None)
               if tpl:
                  f = tpl[1]
                  f <<= tpl[0]
                  f.setName(f"<IF True THEN {str(tpl[0])}>")
               elif i.output_definition.has_grammar():
                  raise KeyError(f"In Switch, the item {i.name} for case {value} was not prepared")

           for i in set(self.all_values()).difference(ok.values()):
               tpl = grammars.get(i.name, None)
               if tpl:
                  f = tpl[1]
                  f.setName(f"<IF False THEN {str(tpl[0])}>")
                  f <<= self.empty
               elif i.output_definition.has_grammar():
                  raise KeyError(f"In Switch, the item {i.name} for case {value} was not prepared")

       if not self.container.force_order:
            return grammar
       return grammar.addParseAction(lambda x: choose(x[0][1]) and x)

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

   def switch_value(self, container):
       return self.fn(container.get(self.item))

   @property
   def condition_name(self):
       return self.item


class SeparatorDefinition(VirtualDefinition):
    """ Basic class for separators """
    def __repr__(self):
        return "<SEPARATOR>"

    def __init__(self, separator_type=None, length=None, write_condition=None):
        super().__init__(template='SEPARATOR')
        self.write_condition = write_condition
        if separator_type is not None:
            if length is not None:
                separator_type = Separator(char=separator_type,length=length)
            self.separator_type = separator_type

    def _create_grammar(self, allow_dangerous=False):
        return pp.Suppress(self.separator_type.grammar())

    def _save_to_file(self, file, value, always=False, name_in_grammar=None, delimiter=''):
        if self.write_condition and not self.write_condition(value):
            return False
        if not always:
            if self.condition and self.condition(value):
                    return False
        if delimiter:
            file.write(delimiter)
        self.separator_type.write(file, None)
        return True
