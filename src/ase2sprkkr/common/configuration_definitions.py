"""
Configuration definitions are classes, that desribes
the syntax of a configuration file, or its parts
(sections or configuration options)

They are able both to parse a file, which results in an
instance of (an instance of :py:class:`ase2sprkkr.common.BaseConfiguration`,
e.g. an :py:class:`Option<ase2sprkkr.common.options.Option>` or
:py:class:`Section<ase2sprkkr.common.configuration_containers.Section>`
), or write such object to a file.
"""

from ..common.grammar_types  import type_from_type, type_from_value, BaseType
from ..common.grammar  import delimitedList, end_of_file, generate_grammar
import pyparsing as pp
from ..common.misc import OrderedDict
from .configuration_containers import Section
from .options import Option
import numpy as np
import inspect
from .misc import cache
import itertools

#:This serves just for dealing with various pyparsing versions
_parse_all_name = 'parse_all' if \
  'parse_all' in inspect.getfullargspec(pp.Or.parseString).args \
  else 'parseAll'

def unique_dict(values):
    """ Create a dictionary from the arguments. However, raise an exception,
    if there is any duplicit key """
    out = dict(values)
    if len(values) != len(out):
       seen = set()
       duplicates = {x for x,y in values if x not in seen and seen.add(x)}
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

   def __init__(self, name, alternative_names=None, is_optional=False, is_hidden=False,
                name_in_grammar=None, help=None, description=None, write_alternative_name=False):
       """
       Parameters
       ----------
        name: str
          Name of the value/section

        alternative_names: str or [str]
          Alternative names that can denotes the value

        is_hidden: boolean
          Hidden values are not offered to a user, usually they are
          set by another object (and so a direct setting of their values
          has no sense)

        name_in_grammar: boolean or None
          If False, there the name of the variable is not printed in the
          configuration file. The variable is recognized by its position.
          If None, the default class value is used

        is_optional: boolean
          If True, this section/value can be missing in the .pot/task file

        help: str
          A short help message for the value/section. It will be the perex for description.

        description: str
           The additional informations for the users.
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
       """ Is it required part of configuration (or can it be ommited)? """
       self.write_alternative_name = write_alternative_name
       self.name_in_grammar = self.__class__.name_in_grammar \
                               if name_in_grammar is None else name_in_grammar
       self.help = help
       """ A short help text describing the content for the users. """
       self.description = description
       """ A longer help text describing the content for the users. """

   def create_object(self, container=None):
       """ Creates Section/Option/.... object (whose properties I define) """
       return self.result_class(self, container)

   def grammar(self):
       """ Generate grammar with the correct settings of pyparsing """
       with generate_grammar():
         return self._grammar()

   def _tuple_with_my_name(self, expr, delimiter=None, has_value=True):
        """ Create the grammar returning tuple (self.name, <result of the expr>) """
        if self.name_in_grammar:
            name = pp.CaselessKeyword(self.name)
            if self.do_not_skip_whitespaces_before_name:
               name.leaveWhitespace()

            if self.alternative_names:
               alt_names = ( pp.CaselessKeyword(n) for n in self.alternative_names )
               if self.do_not_skip_whitespaces_before_name:
                  alt_names = ( i.leaveWhitespace() for i in alt_names )
               name = name ^ pp.Or(alt_names)
               if self.do_not_skip_whitespaces_before_name:
                  name.leaveWhitespace()
            name.setParseAction(lambda x: self.name)
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

class BaseValueDefinition(BaseDefinition):

  result_class = Option

  name_in_grammar = None

  def __init__(self, name, type, default_value=None, alternative_names=None,
               fixed_value=None, required=None, help=None, description=None,
               is_hidden=False, name_in_grammar=None, name_format=None,
               is_optional=None):
    """
    Creates the object

    Parameters
    ----------
    name: str
      Name of the configuration value

    type: BaseType
      Configuration value type

    default: mixed
      Default value for configuration

    alternative_names: str or [str]
      Value can have an alternative name (that alternativelly denotes the value)

    fixed: mixed
      If given, this option is not user given, but with fixed_value value (provided by this parameter)

    required: bool
      Required option can not be set to None (however, a required one
      can be still be optional, if it has a default values).
      If required = None, required = not is_optional and default_value is None

    name_in_grammar: bool or None
      The value in the conf file is prefixed by <name><name_value_delimiter>
      If None, the default type value (type.name_in_grammar) is used

    is_optional: bool or None
      If True, the value can be omited, if the fixed order is required
      None means True if required is False

    is_hidden: bool
      The value is hidden from the user (no container.name access to the value)

    name_format: str or None
      The way how the name is written
    """
    self.type = type_from_type(type)
    if default_value is None and not isinstance(self.type, BaseType):
       self.type = type_from_value(type)
       self.default_value = self.type.convert(type)
    else:
       self.default_value = self.type.convert(default_value) if default_value is not None else None
    assert isinstance(self.type, BaseType), "grammar_type (sprkkr.common.grammar_types.BaseType descendat) required as a value type"

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    if required is None:
       required = not is_optional and default_value is None

    if is_optional is None:
       is_optional = required is False

    super().__init__(
         name = name,
         alternative_names = alternative_names,
         is_optional = is_optional,
         is_hidden = is_hidden,
         name_in_grammar = name_in_grammar,
         help = help,
         description = description
    )

    if self.name_in_grammar is None:
        self.name_in_grammar = self.type.name_in_grammar

    self.fixed_value = self.type.convert(fixed_value) if fixed_value is not None else None
    self.required = self.default_value is not None if required is None else required
    self.help = None
    self.is_hidden = is_hidden
    self.is_optional = is_optional
    self.name_format = name_format

  @property
  def formated_name(self):
    name = next(iter(self.alternative_names)) if self.write_alternative_name else self.name
    if self.name_format:
       return "{:{}}".format(name, self.name_format)
    return name

  def validate(self, value):
    if value is None:
       if self.required:
          raise ValueError(f"The value is required for {self.name}, cannot set it to None")
       return True
    if self.fixed_value is not None and not np.array_equal(self.fixed_value, value):
       raise ValueError(f'The value of {self.name} is required to be {self.fixed_value}, cannot set it to {value}')
    self.type.validate(value, self.name)

  @property
  def value_name_format(self):
    return self.name_format

  @value_name_format.setter
  def value_name_format(self, value):
      self.name_format = value

  def __str__(self):
    out="SPRKKR({} of {})".format(self.name, str(self.type))
    val = self.get_value()
    if val is not None:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar(self):
    body = self.type.grammar(self.name)

    if self.fixed_value is not None:
      def check_fixed(s, loc, x, body=body):
          if self.fixed_value.__class__ is np.ndarray:
             eq = np.array_equal(x[0], self.fixed_value)
          else:
             eq = x[0]==self.fixed_value
          if eq:
             return x
          message="The value of {} is {} and it should be {}".format(self.name, x[0], self.fixed_value)
          raise pp.ParseException(s,loc,message, body)
      body=body.copy().addParseAction(check_fixed)

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

    out = self._tuple_with_my_name(body, has_value=self.type.has_value)
    out.setName(self.name + nbody)
    return out

  def get_value(self, option=None):
     if self.fixed_value is not None:
        return self.fixed_value
     if self.default_value is not None:
        if callable(self.default_value):
           return self.default_value(option)
        return self.default_value
     return None

  def write(self, f, value):
     if self.type.has_value:
        if value is None:
           value = self.get_value()
        if value is None:
           return
        missing, df, np = self.type.missing_value()
        if np.__class__ is value.__class__ and np == value:
           return
        write_value = not ( missing and df == value )
     else:
        value=None
        write_value=True
     f.write(self.prefix)
     if self.name_in_grammar:
        f.write(self.formated_name)
        if write_value:
           f.write(self.name_value_delimiter)
           self.type.write(f, value)
     elif write_value:
        self.type.write(f, value)
     return True


  def copy(self, **kwargs):
     if not '_init_args' in self.__class__.__dict__:
        self.__class__._init_args = inspect.getfullargspec(self.__init__).args[1:]
     default = { k: getattr(self, k) for k in self._init_args }
     default.update(kwargs)
     return self.__class__(**default)

  def remove(self, name):
     del self.section[name]
     return self

def add_excluded_names_condition(element, names):
    """ Add the condition to the element, that
    its value is not any of given names """
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class BaseContainerDefinition(BaseDefinition):
    """ Base class for a definition of the format of a container """

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

    def __init__(self, name, members=[], alternative_names=[], help=None, description=None, is_hidden=False, has_hidden_members=False, is_optional=False, name_in_grammar=None, force_order=None):
       super().__init__(
           name = name,
           alternative_names = alternative_names,
           is_optional = is_optional,
           is_hidden = is_hidden,
           name_in_grammar = name_in_grammar,
           help = help,
           description = description
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

    _init_args = inspect.getfullargspec(__init__).args[3:]

    def copy(self, args=[], items=[], remove=[], defaults={}, **kwargs):
        """ Copy the section with the contained values modified by the arguments."""
        members = OrderedDict( ( (k,i.copy()) for k,i in self._members.items() ) )
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = members[i].type.convert(v)

        default = { k: getattr(self, k) for k in self._init_args }
        default.update(kwargs)
        return self.__class__(self.name, members=members, **default)

    def create_object(self, container=None):
        return self.result_class(self, container)

    def all_member_names(self):
        return itertools.chain(
          (i.name for i in self),
          itertools.chain.from_iterable((
            i.alternative_names for i in self if i.alternative_names
          ))
        )

    def _values_grammar(self, delimiter=None):
       if self.custom_class:
          custom_value = self.custom_member_grammar(self.all_member_names())
       else:
          custom_value = None
       delimiter = delimiter or self.grammar_of_delimiter()

       def grammars():
         """ This function iterates over the items of the container, joining all the without name_in_grammar with the previous ones. """
         it = iter(self._members.values())
         head = next(it)
         curr = head._grammar()

         for i in it:
             if i.name_in_grammar:
               yield head, curr
               head = i
               curr = i._grammar()
             else:
               add = delimiter + i._grammar()
               if i.is_optional:
                  add = pp.Optional(add)
               curr = curr + add
         yield head, curr

       if self.force_order:
           if custom_value:
              cvs = pp.ZeroOrMore(custom_value + delimiter)
           else:
              cvs = None
           first = True
           def not_first(x):
               nonlocal first
               first = False
               return x

           inter = pp.Empty().addCondition(lambda x: first).setName('<is_first>') | \
                   ((delimiter + cvs) if cvs else delimiter)
           inter.setName('<is first>|[<custom section>...]<delimiter>')

           def sequence():
               for head,g in grammars():
                   g.addParseAction(not_first)
                   g = inter + g
                   if head.is_optional:
                      g = pp.Optional(g)
                   yield g

           values  = pp.And([ i for i in sequence()])
           if custom_value:
              if not self._first_section_is_fixed():
                  values = cvs + values
              values += pp.ZeroOrMore(delimiter + custom_value)
       else:
           it = grammars()
           #store the first fixed "chain of sections"
           first = self._first_section_is_fixed() and next(it)[1]
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
              is_ok = self.validate(value[0])
              if is_ok is not True:
                raise pp.ParseException(s, loc, is_ok)
              return value
          values.addParseAction(_validate)

       return values

    def _grammar(self):
       delimiter = self.grammar_of_delimiter()
       values = self._values_grammar(delimiter)
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

    custom_name_characters = pp.alphanums + '_-'
    """ Which characters can appears in an unknown child (value/section) name """

    @classmethod
    def custom_member_grammar(cls, value_names = []):
       """ Grammar for the custom - unknown - child """
       name = pp.Word(cls.custom_name_characters).setParseAction(lambda x: x[0].strip())
       add_excluded_names_condition(name, value_names)
       out = (name + cls.delimited_custom_value_grammar()).setParseAction(lambda x: tuple(x))
       out.setName(cls.custom_value_name)
       return out


    def _first_section_is_fixed(self):
       """ Has/ve the first child(s) in an unordered sequence fixed position? """
       return not self._members.first_item().name_in_grammar

    def parse_file(self, file, return_value_only=True):
       """ Parse the file, return the parsed data as dictionary """
       return self.parse_return(self.grammar().parseFile(file, **{ _parse_all_name: True } ), return_value_only)

    def parse(self, str, whole_string=True, return_value_only=True):
       """ Parse the string, return the parsed data as dictionary """
       return self.parse_return(self.grammar().parseString(str, **{ _parse_all_name: whole_string} ), return_value_only)

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

    async def parse_from_stream(self, stream, up_to, start=None, whole_string=True, return_value_only=True):
        """
        Parse string readed from asyncio stream.
        The stream is readed up to the given delimiter
        """

        result = await stream.readuntil(up_to)
        result = result[:-len(up_to)].decode('utf8')
        if start:
           result = start + result
        return self.parse(result, whole_string)



class BaseSectionDefinition(BaseContainerDefinition):
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


class ConfDefinition(BaseContainerDefinition):
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
           if not isinstance(defs, BaseSectionDefinition):
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

   def read_from_file(self, file, **kwargs):
       """ Read a configuration file and return the parsed Configuration object """
       out = self.result_class(definition = self, **kwargs)
       out.read_from_file(file)
       return out

   def _tuple_with_my_name(self, expr, delimiter=None):
       """ Do not create tuple (name, value) for the root class. """
       return expr

   def parse_return(self, val, return_value_only=True):
        """ There is no name in the parsed results (see how
            ConfDefinition._tuple_with_my_name is redefined)
        """
        val = val[0]
        return val

   def _grammar(self):
       """Returns the grammar to parse the configuration file.

       This method just tweaks the grammar (generated by the common container implementation) to ignore comments,
       so the comments would be ignored just once.
       """

       out=super()._grammar()
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out
