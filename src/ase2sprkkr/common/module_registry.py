"""Reusable discovery registry for values stored in Python modules."""

from __future__ import annotations

import importlib
import pkgutil


class RegistryEntryNotFound(KeyError):
    """A requested registry entry is not provided by the package."""


class ModuleRegistry:
    """Lazily discover and load one named value from package modules.

    Module names are discovered with :mod:`pkgutil` without importing the
    modules. Values are imported and materialized only when requested.
    Subclasses can customize extraction, validation and aliasing without
    duplicating discovery machinery.

    ``attribute`` names the value exported by each registry module. The base
    registry returns that value unchanged; whether an exported callable is a
    factory is domain-specific policy and belongs in a specialized registry.
    """

    def __init__(self, package, attribute=None, *, aliases=None, include=None):
        self.package = package
        self.attribute = attribute
        self.aliases = {
            self.normalize_name(alias): self.normalize_name(target)
            for alias, target in (aliases or {}).items()
        }
        self._include = include or (lambda name: not name.startswith("_"))
        self._module_names = None
        self._modules = {}
        self._values = {}

    @staticmethod
    def normalize_name(name):
        return str(name).upper()

    def _discover_module_names(self):
        if self._module_names is None:
            discovered = {}
            for item in sorted(
                pkgutil.iter_modules(self.package.__path__),
                key=lambda item: item.name,
            ):
                if not self._include(item.name):
                    continue
                name = self.normalize_name(item.name)
                if name in discovered:
                    raise ValueError(
                        f"Duplicate module registry key {name!r} in "
                        f"{self.package.__name__}"
                    )
                discovered[name] = item.name
            self._module_names = discovered
        return self._module_names

    def _load_module(self, name):
        if name in self._modules:
            return self._modules[name]

        module_name = self._discover_module_names().get(name)
        if module_name is None:
            raise RegistryEntryNotFound(name)

        module = importlib.import_module(f".{module_name}", self.package.__name__)
        self._modules[name] = module
        return module

    def value_from_module(self, name, module):
        """Extract the registered value from a loaded module."""
        if self.attribute is None:
            raise NotImplementedError
        try:
            return vars(module)[self.attribute]
        except KeyError as exc:
            raise RegistryEntryNotFound(name) from exc

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
        requested = self.normalize_name(name)
        canonical = self.aliases.get(requested, requested)
        return canonical in self._discover_module_names()

    def names(self, *, include_aliases=True):
        """Return discovered names without importing registry modules."""
        module_names = self._discover_module_names()
        names = list(module_names)
        if include_aliases:
            names.extend(
                alias
                for alias, target in self.aliases.items()
                if alias not in module_names and target in module_names
            )
        return tuple(names)

    def values(self, *, include_aliases=False):
        """Load and return registered values."""
        return tuple(self[name] for name in self.names(include_aliases=include_aliases))

    @property
    def modules(self):
        """Load and return canonical registry modules."""
        for name in self.names(include_aliases=False):
            self._load_module(name)
        return dict(self._modules)
