"""
When the data are validated, there can be need for merge two data sources:
the old section and the newly set data. This classes handles the issue.
"""
from . import options


class SectionAdaptor:
  """ This class wraps a container to behave as a read-only dict with
  some other "addons" It is used during validation of a container.
  """

  def __init__(self, container):
      self.container=container

  def __hasitem__(self, name):
      return self.container__hasitem__(name)

  def __getitem__(self, name):
      return self.container.__getitem__(name)()

  def get(self, name, default=None):
      try:
          return self.container.get(name)
      except Exception:
          return default

  def is_dangerous(self, name):
      return getattr(self.container, name).is_dangerous()

  def __repr__(self):
      return f"<Adaptor for {self.container}>"


class MergeSectionDefinitionAdaptor:
    """ This class returns a read-only dict-like class
    that merge values from a dict (e.g. newly parsed data) and from the
    definition of a section """

    def __init__(self, values, definition):
        self.values = values
        self.definition = definition

    def __hasitem__(self, name):
        return self.values.__hasitem__(name) or \
               self.definition.__hasitem__(name)

    def __getitem__(self, name):
        try:
            return self.values[name]
        except KeyError:
            return self.definition[name].get_value()

    def get(self, name, default=None):
        try:
            return self.values[name]
        except KeyError:
            try:
                  return self.definition[name].get_value()
            except KeyError:
                  return default

    def is_dangerous(self, name):
        if name in self.values:
            return isinstance(self.values, options.DangerousValue)
        return False

    def __repr__(self):
        return f'Section {self.definition.name} with values {self.values}'


class MergeSectionAdaptor:
    """ This class returns a read-only dict-like class
    that merge values from a dict (e.g. newly parsed data) and from the
    a section """

    def __init__(self, values, section):
        self.values = values
        self.section = section

    def __hasitem__(self, name):
        return self.values.__hasitem__(name) or \
               self.section.__hasitem__(name)

    def __getitem__(self, name):
        try:
            return self.values[name]
        except KeyError:
            return self.section[name]()

    def get(self, name, default=None):
        try:
            return self.values[name]
        except KeyError:
            try:
                  return self.section[name]()
            except KeyError:
                  return default

    def is_dangerous(self, name):
        if name in self.values:
            return isinstance(self.values, options.DangerousValue)
        return getattr(self.section, name).is_dangerous()

    def __repr__(self):
        return f'Section {self.section.name} with added {self.values}'
