"""Lazy registries backed by modules discovered with :mod:`pkgutil`."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
import importlib
import pkgutil
from types import ModuleType
from typing import Generic, TypeVar


T = TypeVar("T")


def _identity(value: str) -> str:
    return value


class ModuleRegistry(Mapping[str, T], Generic[T]):
    """Map modules in a package to lazily loaded values.

    Module names are discovered with :func:`pkgutil.iter_modules` without
    importing the entries. A module is imported only when its value is
    requested. The registry deliberately knows nothing about the value exposed
    by a module; ``loader`` defines that domain-specific contract.
    """

    def __init__(
        self,
        package: ModuleType,
        loader: Callable[[ModuleType], T],
        *,
        normalize: Callable[[str], str] = _identity,
        include: Callable[[str], bool] | None = None,
        aliases: Mapping[str, str] | None = None,
    ) -> None:
        path = getattr(package, "__path__", None)
        if path is None:
            raise TypeError(f"{package.__name__!r} is not a package")

        self._package = package
        self._loader = loader
        self._normalize = normalize
        self._include = include or (lambda _name: True)
        self._aliases = {
            normalize(alias): normalize(target)
            for alias, target in (aliases or {}).items()
        }
        self._modules: dict[str, str] | None = None
        self._values: dict[str, T] = {}

    def _discover(self) -> dict[str, str]:
        if self._modules is not None:
            return self._modules

        prefix = self._package.__name__ + "."
        discovered = sorted(
            pkgutil.iter_modules(self._package.__path__, prefix),
            key=lambda item: item.name,
        )
        modules: dict[str, str] = {}
        for module in discovered:
            name = module.name.rsplit(".", 1)[-1]
            if not self._include(name):
                continue
            key = self._normalize(name)
            if key in modules:
                raise ValueError(
                    f"Duplicate module registry key {key!r} in "
                    f"{self._package.__name__}"
                )
            modules[key] = module.name

        self._modules = modules
        return modules

    def _canonical_key(self, key: str) -> str:
        key = self._normalize(key)
        return self._aliases.get(key, key)

    def __getitem__(self, key: str) -> T:
        canonical = self._canonical_key(key)
        if canonical in self._values:
            return self._values[canonical]

        try:
            module_name = self._discover()[canonical]
        except KeyError:
            raise KeyError(key) from None

        module = importlib.import_module(module_name)
        value = self._loader(module)
        self._values[canonical] = value
        return value

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        return self._canonical_key(key) in self._discover()

    def __iter__(self) -> Iterator[str]:
        return iter(self._discover())

    def __len__(self) -> int:
        return len(self._discover())
