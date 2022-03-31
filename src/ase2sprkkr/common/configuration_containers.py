""" In this file the common containers of configuration values are,
either for task or potential.

Configuration containers are classes, that holds configuration values
and other containers, and are able to write them to a configuration file,
and that are results of parsing of a configuration file.
"""

from ..common.misc import OrderedDict
from ..common.grammar_types import mixed
from .options import Option
import pyparsing as pp
from .base_configuration import BaseConfiguration
import itertools
import re
from typing import Union, Any, Dict

class ConfigurationContainer(BaseConfiguration):
  """ A container for configuration (problem-definition) options and/or sections.

  Options in the configuration (problem-definition) files are grouped to
  sections, sections are then grouped in a configuration file object.
  This is a base class for these containers.

  """

  def __init__(self, definition, container=None):
      """ Create the container and its members, according to the definition """
      super().__init__(definition, container)
      self._members = OrderedDict()
      """
      The members of the container, in a form of ``{obj.name : obj}``
      """
      self._interactive_members = OrderedDict()
      """
      Non-hidden members of the containers, accesible via sanitized names.
      I.e. via names with whitespaces and other special characters replaced by underscore.
      These sanitized names are then used as names for "attributes" of this container, to
      make the members accesible via ``<container>.<member>`` notation.
      """
      for v in definition.members():
          self._add(v.create_object(self))

  def members(self):
      """ Members of the container. I.e. the option of the section, or sections
      of the configuration file e.t.c.

      Returns
      -------
      members: dict
      A dictionary of the shape ``{ name : member }``

      """
      return self._members

  def _get_member(self, name):
      """
      Return the member of the container of a given name.
      It search either the members and interactive_members containers
      """
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
      """
      The members of the container are accesible as attributes of the container, too.
      Either using their normal, or ``sanitized`` names.
      """
      try:
        out = self._get_member(name)
      except AttributeError as e:
        raise AttributeError(f"There is no value with name {name} in {self}.\nMaybe, you want to add a custom value using the add method?") \
              from e
      return out

  def __getitem__(self, name):
      """
      The members of the container are accesible using ``container["member name"]`` notation.
      """
      return self._members[name]

  def __dir__(self):
      """
      Expose the interactive_members in the container attribute listing.
      Interactive_members are the non-hidden members identified by their sanitized names.
      """
      return itertools.chain(self._interactive_members.keys(), super().__dir__())

  def __contains__(self, name):
      """ The check for existence of a member with the given name."""
      return name in self._members

  def clear(self, do_not_check_required=False):
      """
      Erase all values (or reset to default) from all options in the container
      (ad subcontainers)

      Parameters
      ----------
      do_not_check_required: bool
      Do not check validity of the values after clearing. If ``False`` (default)
      is passed as this argument, the required option without a default value
      (or a section containing such value) throw an exception, which prevents the
      clearing (neverthenless, previous values in the section will be cleared anyway).
      """
      for i in self._members.values():
          i.clear(do_not_check_required)

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
         return self.as_dict()
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


  def set(self, values:Union[Dict[str,Any],str,None]={}, value=None, *, unknown='find', **kwargs):
      """
      Set the value(s) of parameter(s). Usage:

      > input_parameters.set({'NITER': 5, 'NE': [10]})
      or
      > input_parameters.set(NITER=5, NE=[10])

      Parameters
      ----------

      values:
        Dictionary of values to be set, or the name of the value, if the value is given.

      value:
        Value to be set. Setting this argument require to pass string name to the values argument.

      unkwnown: 'add', 'find' or None
        How to handle unknown (not known by the definition) parameters.
        If 'find', try to find the values in descendant containers.
        If 'add', add unknown values as custom values.
        If None, throw an exception.
        Keyword only argument.

      **kwargs: dict
        The values to be set (an alternative syntax as syntactical sugar)
      """
      if values.__class__ is str:
         values = { values : value }
      elif value is not None:
         raise ValueError("If value argument of Container.set method is given,"
         " the values have to be string name of the value")

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
      """ Iterate over all members of the container """
      yield from self._members.values()

  def as_dict(self):
      """
      Return the content of the container as a dictionary.
      Nested containers will be transformed to dictionaries as well.

      """
      out = OrderedDict()
      for i in self:
          value = i.as_dict()
          if value is not None:
              out[i.name] = value
      return out or None

  def to_string(self, *, validate=False):
      """
      Return the configuration (problem definition) in a string.

      Returns
      -------
      configuration:str
      The configuration, as it should be saved in a configuration/problem definition file.
      """
      from io import StringIO
      s = StringIO()
      self.save_to_file(s, validate=validate)
      return s.getvalue()

  def _find_value(self, name):
      """
      Find a value of a given name in self or in any
      of owned subcontainers.

      Parameters
      ----------
      name: str
      A name of the sought options

      Returns
      -------
      value:typing.Optional[ase2sprkkr.common.options.Option]
      The first option of the given name, if such exists. ``None`` otherwise.
      """
      for i in self:
          if i._definition.is_hidden:
             continue
          out = i._find_value(name)
          if out:
             return out

  @staticmethod
  def _interactive_member_name(name):
      """ Create a sanitized name from a member-name.

      The sanitized names are keys in ``interactive_members`` array, and thus
      the members are accesible by ``<container>.<member>`` notation.
      """
      return re.sub(r'[-\s.]','_',name)

  def _add(self, member):
      name = member.name
      self._members[name] = member
      if not member._definition.is_hidden:
          iname = self._interactive_member_name(name)
          if not iname in self._interactive_members:
              self._interactive_members[iname] = member

  def save_to_file(self, file, *, validate=True):
      """ Save the configuration to a file. The method is implemented in the descendants.

      TODO
      ----
      The implementations in the descendant could be probably merged.
      """

      raise NotImplemented()


class BaseSection(ConfigurationContainer):
  """ A section of SPRKKR configuration - i.e. part of the configuration file. """

  def __setattr__(self, name, value):
      """ Setting the (unknown) attribute of a section sets the value of the member
      with a given name """
      if name[0]=='_':
        super().__setattr__(name, value)
      else:
        val = self._get_member(name)
        val.set(value)

  def has_any_value(self) -> bool:
      """
      Return True if any member of the section has value.

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

class Section(BaseSection):
  """ A standard section of a task or potential (whose content is predefinded by SectionDefinition) """

  @property
  def definition(self):
      """ The definition of the section.

      Returns
      -------
      ase2sprkkr.common.configuration_definitions.BaseContainerDefinition
      The definition of the section. I.e. the object that defines, which configuration values
      are in the section, their default values etc.
      """
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
  """ Base class for data of configuration/problem-definition files

  In addition to container capabilities, it can read/write its data from/to file.
  """

  def save_to_file(self, file, *, validate=True):
      """ Save the configuration to a file in a format readable by SPR-KKR i

      Parameters
      ----------
      file: str or file
        File to read the data from

      validate: bool
        Validate the data in the container first and raise an exception,
        if there is an error (e.g. the the data are not complete)
      """
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
      """ Read data from a file

      Parameters
      ----------
      file: str or file
        File to read the data from

      clear_first: bool
        Clear the container first.
        Otherwise, the data in the sections that are not present in the
        file are preserved.
      """
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
