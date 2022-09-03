""" The classes for storing one configuration value. """
from typing import Union
from ..common.grammar_types import mixed
from .base_configuration import BaseConfiguration

class Option(BaseConfiguration):
  """ Class for one option (a configuration value) of SPRKKR - either to
  be used as a part of InputParameters or Potential configuration.
  Usage:

  >>> from ase2sprkkr.sprkkr.calculator import SPRKKR
  >>> calculator = SPRKKR()
  >>> conf = calculator.input_parameters
  >>> conf.ENERGY.ImE = 5.
  >>> conf.ENERGY.ImE()
  5.0
  >>> conf.ENERGY.ImE.info
  'Configuration value ImE'
  >>> conf.ENERGY.ImE.help()                     # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  Configuration value ImE
  Data description
  ----------------
  ...
  """
  def __init__(self, definition, container=None, value=None):
      """"
      Parameters
      ----------
      definition: ValueDefinition
          The value type of the option and its format (in potential and/or task file)

      container:
          The container, that owns the object

      value: mixed
          The value of the option.
      """
      super().__init__(definition, container)
      self._value = value
      self._hook = None
      self._definition.type.enrich(self)

  def __call__(self):
      if self._value is not None:
         return self._value
      return self.default_value

  @property
  def default_value(self):
      return self._definition.get_value(self)

  def set(self, value, *, unknown=None):
      """
      Set the value of the option.

      Parameters
      ----------
      value: mixed
        The new value of the option.

      unknown: str or None
        A dummy argument to make the method compatibile with
        ase2sprkkr.sprkkr.common.configuration_containers.ConfigurationContainer.set()
      """
      if value is None:
          return self.clear()
      value = self._definition.type.convert(value)
      self._definition.validate(value)
      self._value = value
      if self._hook:
        self._hook(self)

  def get(self):
      """ Return the value of self """
      return self()

  def clear(self, do_not_check_required=False):
      """ Clear the value: set it to None """
      if not self._definition.type.has_value:
         return
      if self._definition.default_value is None and not do_not_check_required and self._definition.required:
         raise ValueError(f'Option {self._get_path()} must have a value')
      self._value = None
      if self._hook:
        self._hook(self)

  def save_to_file(self, file, *, validate=True):
      """ Write the name-value pair to the given file, if the value
      is set. """
      d = self._definition
      if not d.type.has_value:
         return d.write(file, None)
      value = self()
      if value is not None:
        if d.is_expert and d.type.is_the_same_value(value, d.default_value):
           return
        return d.write(file, value)
      elif validate and not self._definition.is_optional:
        name = self._get_root_container()
        raise Exception(f'Value {self._get_path()} is None and it is not an optional value. Therefore, I cannot save the {name}')

  @property
  def name(self):
      return self._definition.name

  def as_dict(self, only_changed: Union[bool,str]='basic'):
      d = self._definition
      if only_changed and (only_changed!='basic' or d.is_expert):
           v,c = self.value_and_changed()
           return v if c else None
      return self()

  def value_and_changed(self):
      """ Return value and whether the value was changed

          Returns
          -------
          value:mixed
            The value of the options

          changed:bool
            Whether the value is the same as the default value or not
      """
      d = self._definition
      default = self.default_value
      if self._value is not None:
         return self._value, not d.type.is_the_same_value(self._value, default)
      return default, False

  def _find_value(self, name):
      if self._definition.name == name:
         return self

  def get_path(self):
      return self._get_path()

  def __repr__(self):
      v = self._value
      if v is None:
         v = self.default_value
         if v:
            o=' (default)'
            v=' = '+str(v)
         else:
            o='out'
            v=''
      else:
          v=' = ' + str(v)
          o=''

      return f"<Option {self._get_path()} of type {self._definition.type} with{o} value{v}>"

class CustomOption(Option):
  """ An user-added option (configuration value). It can be removed from the section. """

  def remove(self):
      """ Remove me from my "parent" section """
      self._container.remove_member(self._definition.name)

  @classmethod
  def factory(cls, value_definition, type = mixed):
      """ Returns factory function for the given value definition """

      def create(name, section):
          definition = value_definition(name, type)
          definition.removable = True
          return cls(definition, section)

      create.grammar_type = type
      return create
