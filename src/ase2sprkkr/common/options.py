from ..common.grammar_types import mixed
from .conf_common import ConfCommon

class Option(ConfCommon):
  """ Class for one option (a configuration value) of SPRKKR - either to be used as a part of Task or Potential configuration. Usage:

  conf.ENERGY.ImE = 5
  conf.ENERGY.ImE()
  > 5
  conf.ENERGY.ImE.__doc__
  > The ImE option means ...
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

  def __call__(self):
      if self._value is not None:
         return self._value
      return self._definition.get_value(self)

  def set(self, value):
      if value is None:
          return self.clear()
      value = self._definition.type.convert(value)
      self._definition.validate(value)
      self._value = value
      if self._hook:
        self._hook(self)

  @property
  def __doc__(self):
      return self._definition.help

  def clear(self, do_not_check_required=False):
      if not self._definition.type.has_value:
         return
      if not do_not_check_required and self._definition.required:
         raise ValueError(f'Option {self._definition.name} have to had a value')
      self._value = None
      if self._hook:
        self._hook(self)

  def save_to_file(self, file):
      if not self._definition.type.has_value:
         return self._definition.write(file, None)
      value = self()
      if value is not None:
        return self._definition.write(file, value)
      elif not self._definition.is_optional:
        name = self._get_root_container().__class__.__name__.lower()
        raise Exception(f'Value {self._get_path()} is None and it is not an optional value. Therefore, I cannot save the {name}')

  @property
  def name(self):
      return self._definition.name

  def to_dict(self, dct):
      value = self()
      if value is not None:
         dct[self.name] = value

  def _find_value(self, name):
      if self._definition.name == name:
         return self

  def get_path(self):
      return self._get_path()

class CustomOption(Option):
  """ An user-added option (configuration value). It can be removed from the section. """

  def remove(self):
      self.section.remove(self._definition.name)

  @classmethod
  def factory(cls, value_definition, type = mixed):
      """ Returns factory function for the given value definition """

      def create(name, section):
          definition = value_definition(name, type)
          definition.removable = True
          return cls(definition, section)

      create.grammar_type = type
      return create