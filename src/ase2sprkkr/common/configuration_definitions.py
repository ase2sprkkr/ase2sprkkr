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

_parse_all_name = 'parse_all' if \
  'parse_all' in inspect.getfullargspec(pp.Or.parseString).args \
  else 'parseAll'

def unique_dict(values):
    out = dict(values)
    if len(values) != len(out):
       seen = set()
       duplicates = {x for x,y in values if x not in seen and seen.add(x)}
       raise pp.ParseException(f"There are non-unique keys: {duplicates}")
    return out

class BaseDefinition:

   """ Redefine in descendants """
   result_class = None

   """ By default, all values,sections.... are named (in the configuration file)"""
   name_in_grammar = True

   def __init__(self, name, alternative_names=None, is_optional=False, is_hidden=False,
                name_in_grammar=None, description=None, help=None, write_alternative_name=False, grammar_of_delimiter=None):
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
       self.write_alternative_name = write_alternative_name
       self.name_in_grammar = self.__class__.name_in_grammar \
                               if name_in_grammar is None else name_in_grammar
       self.help = help
       self.description = description

       if grammar_of_delimiter is not None:
          if not callable(grammar_of_delimiter):
             grammar_of_delimiter = lambda x: grammar_of_delimiter
          self.grammar_of_delimiter = grammar_of_delimiter

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
    if val:
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
    if not names:
       return
    names = set((i.upper() for i in names))
    element.addCondition(lambda x: x[0].upper() not in names)

class BaseContainerDefinition(BaseDefinition):

    """ Force order of its members """
    force_order = False

    """ The (print) format, how the name is written """
    value_name_format = None

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
        """ copy the section with the contained values modified """
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
         """ Iterate over grammar, join all with no name_in_grammar with the previous """
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

    """ A function for validation of just the parsed result (not the user input) """
    validate = None


    @classmethod
    @cache
    def delimited_custom_value_grammar(cls):
        return cls.child_class.grammar_of_delimiter() + cls.custom_value_grammar()

    custom_name_characters = pp.alphanums + '_-'

    @classmethod
    def custom_member_grammar(cls, value_names = []):
       name = pp.Word(cls.custom_name_characters).setParseAction(lambda x: x[0].strip())
       add_excluded_names_condition(name, value_names)
       out = (name + cls.delimited_custom_value_grammar()).setParseAction(lambda x: tuple(x))
       out.setName(cls.custom_value_name)
       return out


    def _first_section_is_fixed(self):
       return not self._members.first_item().name_in_grammar

    def parse_file(self, file, return_value_only=True):
      return self.parse_return(self.grammar().parseFile(file, **{ _parse_all_name: True } ), return_value_only)

    def parse(self, str, whole_string=True, return_value_only=True):
      return self.parse_return(self.grammar().parseString(str, **{ _parse_all_name: whole_string} ), return_value_only)

    def parse_return(self, val, return_value_only=True):
        val = val[0]
        if return_value_only:
           val = val[1]
        return val

    async def parse_from_stream(self, stream, up_to, start=None, whole_string=True, return_value_only=True):
      result = await stream.readuntil(up_to)
      result = result[:-len(up_to)].decode('utf8')
      if start:
         result = start + result
      return self.parse(result, whole_string)



class BaseSectionDefinition(BaseContainerDefinition):
   """ Base class for sections in Pot or InputParameters file """


   """ Is the sectio named, or it is its position given by position?
       In the second case, a custom value is not allowed between
       the section and its predecessor
   """
   result_class = Section

   @property
   def values(self):
       return self._members


   """ Type for a not-specified value (e.g. flag) """
   custom_value_name = 'CUSTOM_VALUE'

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

   name_in_grammar = False
   """ No data are given just by a section name - no Flag equivalent for sections """
   optional_value = None

   @classmethod
   def from_dict(cls, name, defs=None):
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

   def read_from_file(self, file, **kwargs):
       out = self.result_class(definition = self, **kwargs)
       out.read_from_file(file)
       return out

   def _tuple_with_my_name(self, expr, delimiter=None):
       """ Do not create tuple (name, value) for the root class """
       return expr

   def parse_return(self, val, return_value_only=True):
        """ There is no name in the parsed results (see how
            ConfDefinition._tuple_with_my_name is redefined)
        """
        val = val[0]
        return val

   def _grammar(self):
       """Ignore comments"""
       out=super()._grammar()
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out
