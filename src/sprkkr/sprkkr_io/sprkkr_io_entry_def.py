from .sprkkr_io_types  import map_type
import pyparsing as pp
from contextlib import contextmanager

class BaseEntry:
   def grammar(self):
       with generate_grammar():
         return self._grammar()



class EntryValueDef(BaseEntry):

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
    self.type = map_type(type)
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
    suffix = pp.Suppress("=") + self.type.grammar();
    if self.type.is_optional():
      suffix = Optional(suffix).setParseAction( lambda x: x or [ True ])
      nsuffix=''
    else:
      nsuffix='=' + self.type.grammar_name()
    if self.fixed_value is not None:
      suffix.addCondition(lambda x: x[0]==self.fixed_value,
         message="Value {} should be fixed to {}".format(self.name, self.fixed_value)
      )
    out = pp.CaselessKeyword(self.name).setParseAction(lambda x: self.name) + suffix
    out.setParseAction(lambda x: tuple(x))
    out.setName(self.name + nsuffix)
    return out

  def get_value(self):
     if self.fixed_value is not None:
        return self.fixed_value
     if self.default_value is not None:
        return self.default_value
     return None

  def write(f, value):
     if value is None:
        value = self.get_value()
        if value is None:
           return
     f.write("     " + self.name)
     if not self.type.is_optional():
       f.write("=")
       self.type.write(f, value)
     f.write("\n")

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
     del self.sections[name]
     return self

def _dict_from_named_values(args, values=None):
    """auxiliary method that creates dictionary from the arguments"""
    values = values or {}
    for value in args:
       values[value.name] = value
    return values


class EntrySectionDef(BaseEntry):

   def __init__(self, name, *args, values=None):
       self.name = name
       self.values = _dict_from_named_values(args, values=values)

   def copy(self, *args, values=None):
        """ copy the section with the contained values modified """
        values = self.values.copy()
        values.update(_dict_from_named_values(args, values))
        return ConfigurationSection(self.name, values=values)

   def _grammar(self):
        out = pp.CaselessKeyword(self.name)

        anyvalue = pp.MatchFirst((i._grammar() for i in self.values.values()))
        values = pp.OneOrMore(pp.Optional(pp.Suppress(pp.LineEnd())) + anyvalue)
        values.setParseAction(lambda x: dict(x.asList()))

        out = (
            pp.CaselessKeyword(self.name).setParseAction(lambda x: self.name) +
            values
        )
        out.setParseAction(lambda x: tuple(x))
        out.setName(self.name)
        return out

   def __getitem__(self, key):
       return self.values[key]

   def remove(self, name):
       del self.values[name]
       return self


class EntryDef(BaseEntry):

   @staticmethod
   def from_dict(defs):

       def gen(i):
           section = defs[i]
           if not isinstance(defs, EntrySectionDef):
              section = EntrySectionDef(i, *section)
           return section

       return EntryDef(*[ gen(i) for i in defs])

   def __init__(self, *sections, help=None):
       self.help = help
       self.sections = _dict_from_named_values(sections)

   def _grammar(self):
       anysection=pp.MatchFirst(( i._grammar() for i in self.sections.values() ))
       delimited = anysection + \
            ( pp.Suppress(pp.LineEnd() + pp.OneOrMore(pp.LineEnd())) | pp.StringEnd())
       out = pp.OneOrMore(delimited)
       out.setParseAction(lambda x: dict(x.asList()))
       out.ignore("#" + pp.restOfLine + pp.LineEnd())
       return out

   def __getitem__(self, key):
       return self.sections[key]


@contextmanager
def generate_grammar():
    """ Set the pyparsing and then restore to the original state """
    try:
      old = None
      pe = pp.ParserElement
      if hasattr(pe, "DEFAULT_WHITE_CHARS"):
          old = pe.DEFAULT_WHITE_CHARS
      pe.setDefaultWhitespaceChars(' \t')
      yield
    finally:
      if old is not None:
        pe.setDefaultWhitespaceChars(old)
