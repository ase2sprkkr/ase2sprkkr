from ..common.grammar_types  import type_from_type, type_from_value, BaseType, mixed
from ..common.grammar  import delimitedList, end_of_file, generate_grammar
import pyparsing as pp
from ..common.misc import OrderedDict
from .conf_containers import Section
from .options import Option
import numpy as np

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

   def create_object(self, container=None):
       return self.result_class(self, container)

   def grammar(self):
       """ Generate grammar with the correct settings of pyparsing """
       with generate_grammar():
         return self._grammar()

   def _tuple_with_my_name(self, expr, delimiter=None):
        """ Create the grammar returning tuple (self.name, <result of the expr>) """
        if self.name_in_grammar:
            name = pp.CaselessKeyword(self.name).setParseAction(lambda x: self.name)
            if delimiter:
              name += delimiter
        else:
            name = pp.Empty().setParseAction(lambda x: self.name)
        out = name - expr
        out.setParseAction(lambda x: tuple(x))
        return out



class BaseValueDefinition(BaseDefinition):

  name_in_grammar = True
  result_class = Option

  def __init__(self, name, type, default_value=None,
               fixed_value=None, required=False, help=None,
               is_hidden=False, name_in_grammar=True,
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

    fixed: mixed
      If given, this option is not user given, but with fixed_value value (provided by this parameter)

    is_hidden: boolean
      A hidden value is not offered to a user, usually they are
      set by another another object (and so their direct settings
      has no sense)

    name_in_grammar: boolean
      If False, there the name of the variable is not printed in the configuration file. The variable is recognized by its      position.

    is_optional: boolean
      If True, value can be missing in the .pot file

    required: bool
      Is this option required?
    """
    self.name = name
    self.name_in_grammar = name_in_grammar
    self.type = type_from_type(type)
    if default_value is None and not isinstance(self.type, BaseType):
       self.default_value = type
       self.type = type_from_value(type)
    else:
       self.default_value = default_value

    if self.default_value is None and self.type.default_value is not None:
       self.default_value = self.type.default_value

    self.fixed_value = self.type.convert(fixed_value) if fixed_value is not None else None
    self.required = default_value is not None if required is None else required
    self.help = None
    self.is_hidden = is_hidden
    self.is_optional = is_optional

  def validate(self, value):
    if value is None:
       if self.required:
          raise ValueError(f"The value is required for {self.name}, cannot set it to None")
       return True
    if self.fixed_value and not np.array_equal(self.fixed_value, value):
       raise ValueError(f'The value of {self.name} is required to be {self.fixed_value}, cannot set it to {value}')
    self.type.validate(value, self.name)

  @property
  def optional_value(self):
    return self.name, None

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
          if x[0]==self.fixed_value:
             return x
          message="The value of {} is {} and it should be {}".format(self.name, x[0], self.fixed_value)
          raise pp.ParseException(s,loc,message, body)
      body.addParseAction(check_fixed)

    if self.name_in_grammar:
       god = self._grammar_of_delimiter()
       body = god + body

    def fail(s,loc, expr, err):
        err.__class__ = pp.ParseFatalException
        breakpoint()
        raise err

    optional, df, np = self.type.missing_value()
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

  def get_value(self):
     if self.fixed_value is not None:
        return self.fixed_value
     if self.default_value is not None:
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
        f.write(self.name)
        if write_value:
           f.write(self.name_value_delimiter)
           self.type.write(f, value)
     elif write_value:
        self.type.write(f, value)
     return True

  def copy(self, **kwargs):
     default = {
         'default_value' : self.default_value,
         'fixed_value' : fixed_value,
         'required' : required,
         'help' : help
     }
     default.update(kwargs)
     self.__class__(self.name, self.type, **kwargs)

  def remove(self, name):
     del self.section[name]
     return self

def add_excluded_names_condition(element, names):
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class BaseDefinitionContainer(BaseDefinition):

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or OrderedDict()
        for value in args:
           items[value.name] = value
        return items

    def __init__(self, name, members=[], help=None, description=None, is_hidden=False, has_hidden_members=False, required=True):
       self.help = help
       self.required = True
       self.name = name
       self.description = description
       self.is_hidden = is_hidden
       if not isinstance(members, OrderedDict):
          members = self._dict_from_named_values(members)
       self._members = members
       self.has_hidden_members = has_hidden_members

    def __iter__(self):
        return iter(self._members)

    def members(self):
        return self._members.values()

    def names(self):
        return self._members.keys()

    def __getitem__(self, key):
        return self._members[key]

    def __setitem__(self, key, value):
        self._members[key]=value

    def remove(self, name):
        del self._members[name]
        return self

    def copy(self, args=[], items=[], remove=[], defaults={}, **kwargs):
        """ copy the section with the contained values modified """
        members = self._members.copy()
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = v

        k=['help', 'description', 'is_hidden', 'has_hidden_members']
        for i in k:
            if not i in kwargs:
               kwargs[i] = getattr(self, i)
        if len(kwargs) > len(k):
            missing = set(kwargs.keys()) - set(k)
            raise TypeError(f'No {",".join(missing)} keyword arguments for copying sections')

        return self.__class__(self.name, members=members, **kwargs)

    def create_object(self, container=None):
        return self.result_class(self, container)


class BaseSectionDefinition(BaseDefinitionContainer):
   """ Base class for sections in Pot or Task file """


   """ Is the sectio named, or it is its position given by position?
       In the second case, a custom value is not allowed between
       the section and its predecessor
   """
   name_in_grammar = True
   result_class = Section

   """ Force order of its members """
   force_order = False

   def __init__(self, name, members=[], name_in_grammar=True, **kwargs):
       """
       Parameters
       ----------
       name: str
           Name of the section

       members: iterable
           Members of the section.

       name_in_grammar: boolean
           If True, the section (in the file) begins with its name
           If False, there are only member name - member value in the file

       **kwargs
           See BaseDefinitionContainer.__init__
       """
       super().__init__(name, members, *kwargs)
       self.name_in_grammar = name_in_grammar

   @property
   def values(self):
       return self._members

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
           values = pp.MatchFirst((i._grammar() for i in self.members()))
           values |= custom_value
           values = delimitedList(values, delimiter)

       values.setParseAction(lambda x: unique_dict(x.asList()))
       out = self._tuple_with_my_name(values, delimiter)
       out.setName(self.name)
       return out

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

   def _grammar(self):
       sections = pp.MatchFirst(( i._grammar() for i in self.sections.values() ))
       sections |= self.custom_section( self.sections.keys() )
       delimited = sections + (self.__class__._grammar_of_delimiter() | end_of_file)
       out = pp.OneOrMore(delimited)
       out.setParseAction(lambda x: unique_dict(x.asList()))
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out

   @classmethod
   def custom_section(cls, section_names = []):
       name = pp.Word(pp.alphanums + '_')
       add_excluded_names_condition(name, section_names)
       value = cls.section_class._grammar_of_delimiter() + cls._custom_section_value()
       out = (name + value).setParseAction(lambda x: tuple(x))
       out.setName('CUSTOM SECTION')
       return out

   def read_from_file(self, file):
       out = self.result_class(self)
       out.read_from_file(file)
       return out
