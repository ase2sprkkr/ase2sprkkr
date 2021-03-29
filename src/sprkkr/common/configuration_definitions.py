from ..common.grammar_types  import type_from_type, type_from_value, BaseType, mixed
from ..common.grammar  import delimitedList, end_of_file, generate_grammar
import pyparsing as pp
from ..common.misc import OrderedDict
from .conf_containers import Section
from .options import Option
import numpy as np
import inspect

def unique_dict(values):
    out = dict(values)
    if len(values) != len(out):
       seen = set()
       duplicates = [x for x,y in values if x not in seen and not seen.add(x)]
       raise pp.ParseException(f"There are non-unique keys: {duplicates}")
    return out

class BaseDefinition:

   """ Redefine in descendants """
   result_class = None

   """ By default, all values,sections.... are named (in the configuration file)"""
   name_in_grammar = True

   def __init__(self, name, alternative_names=None, is_optional=False, is_hidden=False,
                name_in_grammar=None, description=None, help=None):
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
          If True, this section/value can be missing in the .poti/task file

        description: str
          Short help message

        help: str
          (Long) help message for the value/section
       """
       self.name = name
       if isinstance(alternative_names, str):
          self.alternative_names = [ alternative_names ]
       else:
          self.alternative_names = alternative_names
       self.is_optional = is_optional
       self.is_required = is_hidden
       self.name_in_grammar = self.__class__.name_in_grammar \
                               if name_in_grammar is None else name_in_grammar
       self.help = help
       self.description = description

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
            if self.alternative_names:
               name = name | pp.Or( pp.CaselessKeyword(n) for n in self.alternative_names)
            name.setParseAction(lambda x: self.name)
            if delimiter:
              name += delimiter
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)
        out = name - expr
        if has_value:
            return out.setParseAction(lambda x: tuple(x))
        else:
            return out.suppress()

class BaseValueDefinition(BaseDefinition):

  result_class = Option

  def __init__(self, name, type, default_value=None, alternative_names=None,
               fixed_value=None, required=False, help=None, description=None,
               is_hidden=False, name_in_grammar=None, name_format=None,
               is_optional=False):
    """
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
      Required option can not be set to None (however, required can be still be optional, if it has
      default values)

    name_in_grammar: bool or None
      The value in the conf file is prefixed by <name><name_value_delimiter>
      If None, the default type value (type.name_in_grammar) is used

    is_optional: bool
      If True, the value can be omited, if fixed order is required

    is_hidden: bool
      The value is hidden from the user (no container.name access to the value)

    name_format: str or None
      The way the name is written
    """
    self.type = type_from_type(type)
    if default_value is None and not isinstance(self.type, BaseType):
       self.default_value = type
       self.type = type_from_value(type)
    else:
       self.default_value = default_value

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    if name_in_grammar is None:
        name_in_grammar = self.type.name_in_grammar

    super().__init__(
         name = name,
         alternative_names = alternative_names,
         is_optional = is_optional,
         is_hidden = is_hidden,
         name_in_grammar = name_in_grammar,
         help = help,
         description = description
    )

    self.fixed_value = self.type.convert(fixed_value) if fixed_value is not None else None
    self.required = default_value is not None if required is None else required
    self.help = None
    self.is_hidden = is_hidden
    self.is_optional = is_optional
    self.name_format = name_format

  @property
  def formated_name(self):
    if self.name_format:
       return "{:{}}".format(self.name, self.name_format)
    return self.name

  def validate(self, value):
    if value is None:
       if self.required:
          raise ValueError(f"The value is required for {self.name}, cannot set it to None")
       return True
    if self.fixed_value is not None and not np.array_equal(self.fixed_value, value):
       raise ValueError(f'The value of {self.name} is required to be {self.fixed_value}, cannot set it to {value}')
    self.type.validate(value, self.name)

  @property
  def optional_value(self):
    return self.name, None

  @property
  def value_name_format(self):
    return self.name_format

  @value_name_format.setter
  def value_name_format(self, value):
      self.name_format = value

  def __str__(self):
    out="SPRKKR({} of {})".format(self.name, str(self.type))
    val = self.get_value()
    if val:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar(self):
    body = self.type.grammar()

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
      body.addParseAction(check_fixed)

    if self.name_in_grammar:
       god = self._grammar_of_delimiter()
       body = god + body

    optional, df, _ = self.type.missing_value()
    if optional:
      body = pp.Optional(body).setParseAction( lambda x: x or df )
      nbody=''
    else:
      nbody=self.type.grammar_name()
      if self.name_in_grammar:
         nbody = str(god) + nbody

    out = self._tuple_with_my_name(body)
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

  _init_args = inspect.getfullargspec(__init__).args[1:]

  def copy(self, **kwargs):
     default = { k: getattr(self, k) for k in self._init_args }
     default.update(kwargs)
     self.__class__(args)

  def remove(self, name):
     del self.section[name]
     return self

def add_excluded_names_condition(element, names):
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class BaseDefinitionContainer(BaseDefinition):

    """ Force order of its members """
    force_order = False

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or OrderedDict()
        for value in args:
           items[value.name] = value
        return items

    def __init__(self, name, members=[], alternative_names=[], help=None, description=None, is_hidden=False, has_hidden_members=False, is_optional=False, name_in_grammar=None):
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
       self._members = members
       self.has_hidden_members = has_hidden_members

    def __iter__(self):
        return self._members.values()

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
        """ copy the section with the contained values modified """
        members = self._members.copy()
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = v

        default = { k: getattr(self, k) for k in self._init_args }
        default.update(kwargs)
        return self.__class__(self.name, members=members, **default)

    def create_object(self, container=None):
        return self.result_class(self, container)

    def _grammar(self):
       custom_value = self.custom_value(self._members.keys())
       delimiter = self._grammar_of_delimiter()
       if self.force_order:
           out = []
           cvs = pp.ZeroOrMore(custom_value + delimiter)
           first = True
           def grammar(i):
               nonlocal first
               was_first = first
               grammar = i._grammar()
               if first:
                  first = False
               else:
                  grammar = delimiter + grammar
               if i.is_optional:
                   """ If the first item is optional, the delimiter is included and the next item
                       is considered to be the first (without the delimiter)
                       Limitation: there can not be a section with just one optional item
                   """
                   if was_first:
                      first = True
                      grammar += delimiter
                   grammar = (grammar | pp.Empty().setParseAction(lambda x: i.optional_value))
               if i.name_in_grammar:
                   grammar = cvs + grammar
               return grammar
           values  = pp.And([ grammar(i) for i in self._members.values()])
       else:
           values = pp.MatchFirst([i._grammar() for i in self.members()])
           values |= custom_value
           values = delimitedList(values, delimiter)

       values.setParseAction(lambda x: unique_dict(x.asList()))
       out = self._tuple_with_my_name(values, delimiter)
       out.setName(self.name)
       return out



class BaseSectionDefinition(BaseDefinitionContainer):
   """ Base class for sections in Pot or Task file """


   """ Is the sectio named, or it is its position given by position?
       In the second case, a custom value is not allowed between
       the section and its predecessor
   """
   result_class = Section

   @property
   def values(self):
       return self._members

   @classmethod
   def custom_value(cls, value_names = []):
      name = pp.Word(pp.alphanums + '_')
      add_excluded_names_condition(name, value_names)
      value = cls.value_class._grammar_of_delimiter() + mixed.grammar()
      if cls.optional_value:
        value |= cls.optional_value.grammar()
      return (name + value).setParseAction(lambda x: tuple(x))

   """ Type for not-specified value (e.g. flag) """
   optional_value = None

   def _custom_value():
      """ Type for value given by user """
      return mixed

class ConfDefinition(BaseDefinitionContainer):

   name_in_grammar = False

   @classmethod
   def from_dict(cls, name, defs=None):
       def gen(i):
           section = defs[i]
           if not isinstance(defs, BaseSectionDefinition):
              section = cls.section_class(i, section)
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

   @classmethod
   def custom_value(cls, section_names = []):
       name = pp.Word(pp.alphanums + '_')
       add_excluded_names_condition(name, section_names)
       value = cls.section_class._grammar_of_delimiter() + cls._custom_section_value()
       out = (name + value).setParseAction(lambda x: tuple(x))
       out.setName('CUSTOM SECTION')
       return out

   def read_from_file(self, file):
       out = self.result_class(definition = self)
       out.read_from_file(file)
       return out

   def _tuple_with_my_name(self, expr, delimiter=None):
       """ Do not create tuple (name, value) for the root class """
       return expr

   def _grammar(self):
       """Ignore comments"""
       out=super()._grammar()
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out
