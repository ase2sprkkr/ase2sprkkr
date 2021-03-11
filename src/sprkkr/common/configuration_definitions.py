from ..common.grammar_types  import type_from_type, type_from_value, BaseType, mixed
from ..common.grammar  import BaseGrammar, delimitedList
import pyparsing as pp
from ..common.misc import OrderedDict

def unique_dict(values):
    out = dict(values)
    if len(values) != len(out):
       seen = set()
       duplicates = [x for x,y in values if x not in seen and not seen.add(x)]
       raise pp.ParseException(f"There are unique KEYS: {duplicates}")
    return out


class BaseValueDefinition(BaseGrammar):

  name_in_grammar = True

  def __init__(self, name, type, default_value=None,
               fixed_value=None, required=False, help=None):
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

    required: bool
    Is this option required?
    """
    self.name = name
    self.type = type_from_type(type)
    if default_value is None and not isinstance(self.type, BaseType):
       self.default_value = type
       self.type = type_from_value(type)
    else:
       self.default_value = default_value

    self.fixed_value = fixed_value
    self.required = required
    self.help = None

  def __str__(self):
    out="SPRKKR({} of {})".format(self.name, str(self.type))
    val = self.get_value()
    if val:
      out+= "={}".format(val)
    return out

  def __repr__(self):
    return str(self)

  def _grammar(self):
    suffix = self.type.grammar()
    if self.name_in_grammar:
       god = self._grammar_of_delimiter()
       suffix = god + suffix

    optional, df, np = self.type.missing_value()
    if optional:
      suffix = pp.Optional(suffix).setParseAction( lambda x: x or df )
      nsuffix=''
    else:
      nsuffix=str(god) + self.type.grammar_name()
    if self.fixed_value is not None:
      suffix.addCondition(lambda x: x[0]==self.fixed_value,
         message="Value {} should be fixed to {}".format(self.name, self.fixed_value)
      )

    out = self._tuple_with_my_name(suffix)
    out.setName(self.name + nsuffix)
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

class BaseDefinitionContainer(BaseGrammar):

    @staticmethod
    def _dict_from_named_values(args, items=None):
        """auxiliary method that creates dictionary from the arguments"""
        items = items or OrderedDict()
        for value in args:
           items[value.name] = value
        return items

    def __init__(self, name, members=[], help=None, description=None):
       self.help = help
       self.name = name
       self.description = description
       if not isinstance(members, OrderedDict):
          members = self._dict_from_named_values(members)
       self._members = members

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

    def copy(self, args=[], items=[], remove=[], defaults={}, help=None, description=None):
        """ copy the section with the contained values modified """
        members = self._members.copy()
        for i in remove:
            del members[i]
        members.update(self._dict_from_named_values(args, items))
        for i,v in defaults.items():
            members[i].default_value = v
        if help is None: help = self.help
        if description is None: description = self.description
        return self.__class__(self.name, members=members, help=help, description=description)




class BaseSectionDefinition(BaseDefinitionContainer):
   """ Base class for sections in Pot or Task file """


   """ Is the sectio named, or it is its position given by position?
       In the second case, a custom value is not allowed between
       the section and its predecessor
   """
   name_in_grammar = True
   """ Force order of its members """
   force_order = False

   @property
   def values(self):
        return self._members

   def _grammar(self):
        out = pp.CaselessKeyword(self.name)
        custom_value = self.custom_value(self._members.keys())
        delimiter = self._grammar_of_delimiter()
        if self.force_order:
           out = []
           cvs = pp.ZeroOrMore(custom_value + delimiter)
           def generator():
               for i in self._members.values:
                   if i.name_in_grammar:
                       yield cvs
                   yield i._grammar()
               yield cvs
           out = pp.An,d(*generator())
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
       delimited = sections + (self.__class__._grammar_of_delimiter() | pp.StringEnd())
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
