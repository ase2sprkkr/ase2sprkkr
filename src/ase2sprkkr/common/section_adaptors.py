"""
When the data are validated, there can be need for merge two data sources:
the old section and the newly set data. This classes handles the issue.
"""

from . import options


class SectionAdaptor:
    """This class wraps a container to behave as a read-only dict with
    some other "addons" It is used during validation of a container.
    """

    def __init__(self, container):
        self.container = container

    def __contains__(self, name):
        return name in self.container

    def __getitem__(self, name):
        return self.container.__getitem__(name)()

    def get(self, name, default=None):
        try:
            return self.container.get(name)
        except Exception:
            return default

    def is_dangerous(self, name):
        return self.container[name].is_dangerous()

    def __repr__(self):
        return f"<Adaptor for {self.container}>"


class MergeSectionDefinitionAdaptor:
    """Read-only view merging parsed data with configuration definitions.

    Nested container values are exposed as adaptors sharing the same ``root``,
    which allows deferred conditions to inspect sibling sections.
    """

    def __init__(self, values, definition, root=None, parent=None):
        self.values = values
        self.definition = definition
        self.parent = parent
        self.root = root or self

    def __contains__(self, name):
        return name in self.values or name in self.definition

    def __getitem__(self, name):
        try:
            definition = self.definition[name]
        except KeyError:
            return self.values[name]
        try:
            value = self.values[name]
        except KeyError:
            if hasattr(definition, "_members"):
                value = {}
            else:
                return definition.get_value()

        if hasattr(definition, "_members"):
            if isinstance(value, list):
                return [self.__class__(item, definition, self.root, self) for item in value]
            return self.__class__(value, definition, self.root, self)
        return value

    def get(self, name, default=None):
        try:
            return self[name]
        except (KeyError, TypeError):
            return default

    def was_parsed(self, name):
        return name in self.values or name in getattr(self.values, "checks", ())

    def is_dangerous(self, name):
        if name in self.values:
            return isinstance(self.values[name], options.DangerousValue)
        return False

    def __repr__(self):
        return f"Section {self.definition.name} with values {self.values}"


class MergeSectionAdaptor:
    """This class returns a read-only dict-like class
    that merge values from a dict (e.g. newly parsed data) and from the
    a section"""

    def __init__(self, values, section):
        self.values = values
        self.section = section

    def __contains__(self, name):
        return name in self.values or name in self.section

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
            return isinstance(self.values[name], options.DangerousValue)
        return self.section[name].is_dangerous()

    def __repr__(self):
        return f"Section {self.section.name} with added {self.values}"
