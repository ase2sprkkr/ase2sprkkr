from ..common.grammar_types import mixed

class Option:
  """ Class for one option (a configuration value) of SPRKKR - either to be used as a part of Task or Potential con figuration. Usage:

  conf.ENERGY.ImE = 5
  conf.ENERGY.ImE()
  > 5
  conf.ENERGY.ImE.__doc__
  > The ImE option means ...
  """
  def __init__(self, definition, value=None):
      """"
      Parameters
      ----------
      definition: ValueDefinition
          The value type of the option and its format (in potential and/or task file)

      value: mixed
          The value of the option.
      """

      self._definition = definition
      self._value = value

  def __call__(self):
      if self._value is not None:
         return self._value
      return self._definition.get_value()

  def set(self, value):
      self._definition.type.validate(value)
      self._value = value

  @property
  def __doc__(self):
      return self._definition.help

  def clear(self):
      self._value = None

  def save_to_file(self, file):
      value = self()
      if value is not None:
        return self._definition.write(file, value)

  @property
  def name(self):
      return self._definition.name

  def to_dict(self, dct):
      value = self()
      if value is not None:
         dct[self.name] = value

class CustomOption(Option):
  """ An user-added option (configuration value). It can be removed from the section. """

  def __init__(self, section, definition):
      super().__init__(definition)

  def remove(self):
      self.section.remove(self._definition.name)

  @classmethod
  def factory(cls, value_definition, type = mixed):
      """ Returns factory function for the given value definition """

      def create(section, name):
          definition = value_definition(name, type)
          definition.removable = True
          return cls(section, definition)
      return create
