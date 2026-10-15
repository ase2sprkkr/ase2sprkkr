"""Reusable discovery registry for values stored in Python modules."""

from __future__ import annotations

import importlib
import pkgutil
from types import MappingProxyType


class RegistryEntryNotFound(KeyError):
    """A requested registry key has no discoverable module."""


class RegistryDefinitionError(RuntimeError):
    """A discovered registry module or registry configuration is malformed."""


class ModuleRegistry:
    """Lazily discover and load one named value from package modules.

    Module names are discovered with :mod:`pkgutil` without importing the
    modules. Values are imported and materialized only when requested.
    Subclasses can customize extraction, validation and aliasing without
    duplicating discovery machinery.

    ``attribute`` names the value exported by each registry module. The base
    registry returns that value unchanged; whether an exported callable is a
    factory is domain-specific policy and belongs in a specialized registry.

    Aliases must point directly to canonical module names. They may not shadow
    modules, form chains, or differ only by normalization. This keeps
    ``names()``, membership tests and item lookup consistent.
    """

    def __init__(self, package, attribute=None, *, aliases=None, include=None):
        self.package = package
        self.attribute = attribute
        self._include = include or (lambda name: not name.startswith("_"))
        self.aliases = MappingProxyType(self._normalize_aliases(aliases or {}))
        self._module_names = None
        self._modules = {}
        self._values = {}

    @staticmethod
    def normalize_name(name):
        return str(name).upper()

    def _normalize_aliases(self, aliases):
        normalized = {}
        for alias, target in aliases.items():
            alias = self.normalize_name(alias)
            target = self.normalize_name(target)
            if alias in normalized:
                raise RegistryDefinitionError(
                    f"Duplicate normalized registry alias {alias!r}"
                )
            normalized[alias] = target
        return normalized

    def _validate_aliases(self, module_names):
        for alias, target in self.aliases.items():
            if alias in module_names:
                raise RegistryDefinitionError(
                    f"Registry alias {alias!r} shadows a module in "
                    f"{self.package.__name__}"
                )
            if target not in module_names:
                raise RegistryDefinitionError(
                    f"Registry alias {alias!r} targets {target!r}, which is not "
                    f"a canonical module in {self.package.__name__}"
                )

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
                    raise RegistryDefinitionError(
                        f"Duplicate module registry key {name!r} in "
                        f"{self.package.__name__}"
                    )
                discovered[name] = item.name
            self._validate_aliases(discovered)
            self._module_names = MappingProxyType(discovered)
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
            raise RegistryDefinitionError(
                f"Registry module {module.__name__!r} does not define required "
                f"export {self.attribute!r}"
            ) from exc

    def prepare_value(self, name, module, value):
        """Validate or transform a freshly loaded value."""
        return value

    def alias_value(self, alias, value):
        """Return the value exposed under *alias*."""
        return value

    def _load_value(self, name):
        if name in self._values:
            return self._values[name]
        module = self._load_module(name)
        value = self.value_from_module(name, module)
        value = self.prepare_value(name, module, value)
        self._values[name] = value
        return value

    def __getitem__(self, name):
        requested = self.normalize_name(name)
        if requested in self._values:
            return self._values[requested]

        module_names = self._discover_module_names()
        canonical = self.aliases.get(requested, requested)
        if canonical not in module_names:
            raise RegistryEntryNotFound(requested)

        value = self._load_value(canonical)
        if requested != canonical:
            value = self.alias_value(requested, value)
            self._values[requested] = value
        return value

    def __contains__(self, name):
        """Return whether *name* is discoverable without importing its module."""
        requested = self.normalize_name(name)
        module_names = self._discover_module_names()
        return requested in module_names or requested in self.aliases

    def names(self, *, include_aliases=True):
        """Return discovered names without importing registry modules."""
        module_names = self._discover_module_names()
        names = list(module_names)
        if include_aliases:
            names.extend(self.aliases)
        return tuple(names)

    def values(self, *, include_aliases=False):
        """Load and return registered values."""
        return tuple(
            self[name]
            for name in self.names(include_aliases=include_aliases)
        )

    @property
    def modules(self):
        """Import and return canonical registry modules without materializing values."""
        for name in self.names(include_aliases=False):
            self._load_module(name)
        return dict(self._modules)
