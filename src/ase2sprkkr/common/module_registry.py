"""Reusable discovery registry for definitions stored in Python modules."""

from __future__ import annotations

import importlib
import pkgutil


class ModuleRegistry:
    """Lazily discover and load one named value from package modules.

    Modules are discovered with :mod:`pkgutil`; values themselves are loaded
    only when requested. Subclasses can customize extraction, validation and
    aliasing without duplicating discovery machinery.
    """

    def __init__(self, package, attribute=None, *, aliases=None):
        self.package = package
        self.attribute = attribute
        self.aliases = {
            self.normalize_name(alias): self.normalize_name(target)
            for alias, target in (aliases or {}).items()
        }
        self._module_names = None
        self._modules = {}
        self._values = {}
        self._missing = set()

    @staticmethod
    def normalize_name(name):
        return str(name).upper()

    def _discover_module_names(self):
        if self._module_names is None:
            self._module_names = {
                item.name.upper(): item.name
                for item in pkgutil.iter_modules(self.package.__path__)
                if not item.name.startswith("_")
            }
        return self._module_names

    def _load_module(self, name):
        if name in self._missing:
            raise KeyError(name)
        if name in self._modules:
            return self._modules[name]

        module_name = self._discover_module_names().get(name)
        if module_name is None:
            self._missing.add(name)
            raise KeyError(name)

        module = importlib.import_module(f".{module_name}", self.package.__name__)
        if self.attribute is not None and not hasattr(module, self.attribute):
            self._missing.add(name)
            raise KeyError(name)
        self._modules[name] = module
        return module

    def value_from_module(self, name, module):
        """Extract the registered value from a loaded module."""
        if self.attribute is None:
            raise NotImplementedError
        value = getattr(module, self.attribute)
        return value() if callable(value) else value

    def prepare_value(self, name, module, value):
        """Validate or transform a freshly loaded value."""
        return value

    def alias_value(self, alias, value):
        """Return the value exposed under *alias*."""
        return value

    def _load_value(self, name):
        module = self._load_module(name)
        value = self.value_from_module(name, module)
        value = self.prepare_value(name, module, value)
        self._values[name] = value
        return value

    def __getitem__(self, name):
        requested = self.normalize_name(name)
        if requested in self._values:
            return self._values[requested]

        canonical = self.aliases.get(requested, requested)
        if canonical in self._values:
            value = self._values[canonical]
        else:
            value = self._load_value(canonical)

        if requested != canonical:
            value = self.alias_value(requested, value)
            self._values[requested] = value
        return value

    def __contains__(self, name):
        try:
            self[name]
        except KeyError:
            return False
        return True

    def names(self, *, include_aliases=True):
        names = []
        for name in self._discover_module_names():
            try:
                self[name]
            except KeyError:
                continue
            names.append(name)
        if include_aliases:
            names.extend(
                alias
                for alias, target in self.aliases.items()
                if target in names
            )
        return tuple(names)

    def values(self, *, include_aliases=False):
        return tuple(self[name] for name in self.names(include_aliases=include_aliases))

    @property
    def modules(self):
        """Loaded registry modules, keyed by normalized name."""
        for name in self.names(include_aliases=False):
            self._load_module(name)
        return dict(self._modules)
