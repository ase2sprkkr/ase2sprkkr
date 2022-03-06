""" In this file the common containers of configuration values are,
either for task or potential """

from ..common.misc import OrderedDict
from ..common.grammar_types import mixed
from .options import Option
import pyparsing as pp
from .base_configuration import BaseConfiguration
import itertools
import re

class ConfigurationContainer(BaseConfiguration):
  """ Custom task section. Section created by user with no definition """

  def __init__(self, definition, container=None):
      super().__init__(definition, container)
      self._members = OrderedDict()
      #to acces the member via container.member notation, the forbidden
      #characters in the name are replaced by _
      self._interactive_members = OrderedDict()
      for v in definition.members():
          self._add(v.create_object(self))

  def members(self):
      return self._members

  def _get_member(self, name):
      if name in self._members:
          out = self._members[name]
      elif name in self._interactive_members:
          out = self._interactive_members[name]
      else:
          raise AttributeError(f'No {name} member of {self._definition}')
      if out._definition.is_hidden:
          raise AttributeError(f'Member {name} of {self._definition} is not directly accessible')
      return out

  def __getattr__(self, name):
      try:
        out = self._get_member(name)
      except AttributeError as e:
        raise AttributeError(f"There is no value with name {name} in {self}.\nMaybe, you want to add a custom value using the add method?") \
              from e
      return out

  def __getitem__(self, name):
    return self._members[name]

  def __iter__(self):
      yield from self._members.values()

  def __dir__(self):
      return itertools.chain(self._interactive_members.keys(), super().__dir__())

  def __contains__(self, name):
      return name in self._members

  def clear(self, do_not_check_required=False):
      for i in self._members.values():
          i.clear(do_not_check_required)

  @property
  def name(self):
      return self._definition.name

  def get(self, name=None, unknown='find'):
      """
      Get the value, either of self or of a child of a given name.

      Parameters
      ----------
      name: None or str
        If None, return contained values as a dictionary.
        Otherwise, return the value of the member with the given name.

      unknown: str or None
        If unknown == 'find' and there is no member with a given name,
        try to find the first such in descendant conainers.

      Return
      ------
      value: mixed
      """

      if name is None:
         return self.to_dict()
      if '.' in name:
         section, name = name.split('.')
         return self._members[section].get(name)
      if name in self._members:
         val = self._members[name]
      elif unknown=='find':
         val = self._find_value(name)
      else:
         val = None
      if not val:
         raise ValueError(f"No {name} member of {self}")
      return val.get()


  def set(self, values={}, *, unknown='find', **kwargs):
      """
      Set the value(s) of parameter(s). Usage:

      > input_parameters.set({'NITER': 5, 'NE': [10]})
      or
      > input_parameters.set(NITER=5, NE=[10])

      Parameters
      ----------

      options: values or None
        Dictionary of values to be set.

      unkwnown: 'add', 'find' or None
        How to handle unknown (not known by the definition) parameters.
        If 'find', try to find the values in descendant containers.
        If 'add', add unknown values as custom values.
        If None, throw an exception.
        Keyword only argument.

      **kwargs: dict
        The values to be set (an alternative syntax as syntactical sugar)
      """

      def set_value(name, value):
        if '.' in name:
          section, name = name.split('.')
          self._members[section].set({name:value})
        if name not in self._members:
           if unknown == 'find':
              option = self._find_value(name)
              if option:
                 option.set(value)
                 return
           if not unknown == 'add':
              raise KeyError("No option with name {} in {}".format(name, str(self)))
              return
           self.add(name, value)
        else:
           self._members[name].set(value, unknown=unknown)

      if values:
        for i,v in values.items():
           set_value(i,v)
      if kwargs:
        for i,v in kwargs.items():
          set_value(i,v)

  def add(self, name:str, value=None):
      """
      Add custom value to the container

      Parameters
      ----------
      name: str
        Name of the added value

      value: value
        Value of the added value
      """
      if not getattr(self._definition, 'custom_class', False):
         raise TypeError(f'Can not add custom members to a configuration class {self._definition}')
      if name in self._members:
         raise TypeError(f'Section member {name} is already in the section {self._definition}')
      cc = self._definition.custom_class
      self._add(cc(name, self))
      if value is not None:
          self._members[name].set(value, unknown='add')


  def remove_member(self, name:str):
      """
      Remove a (previously added) custom value from the container
      """
      cclass = getattr(self._definition, 'custom_class', False)
      if not cclass:
         raise TypeError("Can not remove items of {}".format(name))
      if not getattr(self._members[name], 'remove'):
         raise KeyError("No custom member with name {} to remove".format(name))
      member = self._members[name]
      del self._members[name]
      iname = self._interactive_member_name(name)
      if iname in self._interactive_members and \
          self._interactive_members[iname] == member:
            del self._interactive_members[iname]

  def __iter__(self):
      yield from self._members.values()

  def to_dict(self, dct=None):
      """
      Return the content of the container as a dictionary.
      Nested containers will be transformed to dictionaries as well.
      """
      out = OrderedDict()
      for i in self:
          i.to_dict(out)
      if dct is not None and out:
          dct[self.name] = out
      return out

  def to_string(self, *, validate=False):
      from io import StringIO
      s = StringIO()
      self.save_to_file(s, validate=validate)
      return s.getvalue()

  def _find_value(self, name):
      for i in self:
          if i._definition.is_hidden:
             continue
          out = i._find_value(name)
          if out:
             return out

  @staticmethod
  def _interactive_member_name(name):
      return re.sub(r'[ -]','_',name)

  def _add(self, member):
      name = member.name
      self._members[name] = member
      if not member._definition.is_hidden:
          iname = self._interactive_member_name(name)
          if not iname in self._interactive_members:
              self._interactive_members[iname] = member

  def save_to_file(self, file, *, validate=True):
      raise NotImplemented()


class BaseSection(ConfigurationContainer):
  """ A section of SPRKKR configuration  """

  def __setattr__(self, name, value):
      if name[0]=='_':
        super().__setattr__(name, value)
      else:
        val = self._get_member(name)
        val.set(value)

  def has_any_value(self) -> bool:
      """
      Return
      ------
        has_any_value: bool
            True, if no value in the container is set, False otherwise
      """
      for i in self:
        if i() is not None:
           return True
      return False

  def save_to_file(self, file, *, validate=True):
      """ Save the content of the container to the file (according to the definition)

      Parameters
      ----------
      file: file
        File object (open for writing), where the data should be written
      """
      if not self.has_any_value():
         if validate and not self._definition.is_optional:
            raise ValueError(f"Non-optional section {self._definition.name} has no value to save")
         return
      if self._definition.name_in_grammar:
         file.write(self._definition.name)
         file.write('\n')
      for o in self:
          if o.save_to_file(file, validate=validate):
             file.write(self._definition.delimiter)
      file.flush()

  @property
  def seciton_name(self):
      """
      Name of the section.

      Returns
      -------
        name: str
          The name of the section (according to the definition of the section)
      """

      return self._definition.name


class Section(BaseSection):
  """ A standard section of a task or potential (whose content is predefinded by SectionDefinition) """

  @property
  def definition(self):
      return self._definition


class CustomSection(BaseSection):
  """ Custom task section. Section created by user with no definition """

  def remove(self):
      """ Remove the custom section from the parent container """
      self._container.remove(self.name)

  @classmethod
  def factory(cls, definition_type):
      """ Create a factory for custom values.

      Parameters
      ----------
      definition_type: ase2sprkkr.common.configuration_definitions.BaseDefinition
        Type (definitions) of the custom values created by
        the resulting function

      Return
      ------
      factory: callable
        Factory function of the signature (name: str, container: ase2sprkkr.common.configuration_containers.ConfigurationContainer)
        that created a custom value or section of the given definition

      """
      def create(name, container):
          definition = definition_type(name)
          definition.removable = True
          return cls(definition, container)
      return create


class RootConfigurationContainer(ConfigurationContainer):

  def save_to_file(self, file, *, validate=True):
      """ Save the configuration to a file in a format readable by SPR-KKR """
      if not hasattr(file, 'write'):
         with open(file, "w") as file:
           return self.save_to_file(file, validate=validate)

      it = iter(self)
      i = next(it)
      if i:
        i.save_to_file(file, validate=validate)
        for i in it:
          file.write(self._definition.delimiter)
          i.save_to_file(file, validate=validate)
      file.flush()

  def read_from_file(self, file, clear_first=True):
      values = self._definition.parse_file(file)
      #except Exception as e:
      #   print(e)
      #   breakpoint()
      #   print(e)
      if clear_first:
         self.clear(True)
      self.set(values, unknown='add')

  def find(self, name):
      """ Find a configuration value of a given name in the owned sections """
      return self._find_value(name)
