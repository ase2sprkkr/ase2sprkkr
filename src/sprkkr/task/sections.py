from collections import OrderedDict
from option_types import Mixed

class Option:
  """ Class for one option of SPRKKR. Usage:

  conf.ENERGY.ImE = 5
  conf.ENERGY.ImE()
  > 5
  conf.ENERGY.ImE.__doc__
  > The ImE option means ...
  """
  def __init__(self, definition, value=None):
      self._definition = definition
      self._value = value

  def __call__(self):
      if self._value is not None:
         return self._value
      return self._definition.get_value()

  def set(self, value):
      self._definition.validate(value)
      self._value = value

  @property
  def __doc__(self):
      return self._definition.help

  def clear(self):
      self._value = None

  def save_to_file(self, file):
      value = self()
      if value is not None:
        self._definition.write(file, value)

class BaseSection:
  """ A section of SPRKKR configuration  """

  def __init__(self, definition):
      self._options = OrderedDict()

  def __getattr__(self, name):
      return self._options[name]

  def __getitem__(self, name):
      return self._options[name]

  def __setattr__(self, name, value):
      if name[0]=='_':
        super().__setattr__(name, value)
      else:
        self._options[name].set(value)

  def set(self, values, allow_add=True):
      for i in values:
          if not i in values:
             if not allow_add:
                raise KeyError("No option with name {} in section {}".format(i, self.section_name))
             self.add(i, values[i])
          else:
             self._options[i] = values[i]

  def clear(self):
      for i in self._options:
        i.clear()

  def __iter__(self):
      yield from self._options.values()

  def has_any_value(self):
      for i in self:
        if i() is not None:
           return True
      return False

  def save_to_file(self, file):
      if not self.has_any_value():
         return
      file.write(self.section_name)
      file.write('\n')
      for o in self:
          o.save_to_file(file)
      file.write('\n')

  def add(self, name, value):
      self._options[name] = Option(name, Mixed, value = value)

  def __contains__(self, name):
      return name in self._options

class Section:
  """ A predefined section of SPRKKR task configuration """

  def __init__(self, definition):
      super().__init__(self)
      self._definition = definition
      for d in definition.options:
          self._options[d] = Option(definition.options[d])

  @property
  def definition(self):
      return self._definition

  @property
  def seciton_name(self):
      return self.definition.name


class CustomSection(Section):
  """ Custom task section. Section created by user with no definition """
  def __init__(self, task, name):
      super().__init__(self)
      self.section_name = name
      self._task = task

  def remove_section(self):
      self._task.remove_section(self.section_name)

