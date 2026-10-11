"""Execution profiles for configuring how SPR-KKR calculations are run.

The collection deliberately distinguishes persistent ase2sprkkr profiles from
virtual profiles backed by other configuration sources.  Virtual profiles are
visible through the mapping API, but are never serialized as ase2sprkkr
configuration.
"""

from __future__ import annotations

from collections.abc import Iterator

from .common.repeated_configuration_containers import RepeatedConfigurationContainer
from .common.warnings import DataValidityError


class ReadOnlyProfileValue:
    """Read-only option-like view used by externally backed profiles."""

    def __init__(self, getter, name: str):
        self._getter = getter
        self.name = name

    def __call__(self):
        return self._getter()

    get = __call__

    def set(self, value, **kwargs):
        raise TypeError(f"Profile value {self.name!r} is read-only")

    def __setitem__(self, key, value):
        raise TypeError(f"Profile value {self.name!r} is read-only")

    def __repr__(self):
        return repr(self())


class LegacyProfile:
    """Read-only view of the deprecated global execution options."""

    name = "legacy"
    source = "legacy"
    read_only = True

    def __init__(self, config):
        self.mpi = ReadOnlyProfileValue(lambda: config.running.mpi(), "mpi")
        self.executable_dir = ReadOnlyProfileValue(
            lambda: config.executables.dir(), "executable_dir"
        )
        self.executable_suffix = ReadOnlyProfileValue(
            lambda: config.executables.suffix(), "executable_suffix"
        )

    def as_dict(self, *args, **kwargs):
        return {
            "mpi": self.mpi(),
            "executable_dir": self.executable_dir(),
            "executable_suffix": self.executable_suffix(),
        }

    get = as_dict

    def set(self, *args, **kwargs):
        raise TypeError(
            "The 'legacy' profile is read-only; change the deprecated "
            "running.mpi/executables options or create a named profile instead"
        )


class ProfilesConfigurationContainer(RepeatedConfigurationContainer):
    """Repeated profile section with non-persistent virtual profiles."""

    reserved_names = frozenset(("legacy", "ase"))

    def _legacy_available(self) -> bool:
        config = self._container
        # Presence, rather than difference from the default, matters here: an
        # explicitly configured value equal to its default still is legacy
        # configuration.  ``None`` is represented by clearing an option.
        return any(
            getattr(option, "_value", None) is not None
            for option in (
                config.running.mpi,
                config.executables.dir,
                config.executables.suffix,
            )
        )

    def _virtual(self, name):
        if name == "legacy" and self._legacy_available():
            return LegacyProfile(self._container)
        raise KeyError(name)

    def add(self, id):
        if id in self.reserved_names:
            raise ValueError(f"{id!r} is a reserved profile name")
        return super().add(id)

    def __getitem__(self, name):
        try:
            return self._values[name]
        except KeyError:
            return self._virtual(name)

    def __contains__(self, name):
        return name in self._values or (
            name == "legacy" and self._legacy_available()
        )

    def _virtual_names(self) -> Iterator[str]:
        if self._legacy_available():
            yield "legacy"

    def __iter__(self):
        yield from self._values
        yield from self._virtual_names()

    def __len__(self):
        return len(self._values) + sum(1 for _ in self._virtual_names())

    def items(self):
        yield from self._values.items()
        for name in self._virtual_names():
            yield name, self._virtual(name)

    def values(self):
        yield from self._values.values()
        for name in self._virtual_names():
            yield self._virtual(name)

    # The methods below intentionally use only _values.  Virtual profiles are
    # views of external/legacy state and must never leak into persistence.
    def _as_dict(self, only_changed="basic", generated=False, copy=False):
        out = {}
        for key, value in self._values.items():
            item = value.as_dict(only_changed, generated, copy)
            if item is not None:
                out[key] = item
        return out or None

    def is_changed(self):
        return any(value.is_changed() for value in self._values.values())

    def _save_to_file(self, file, always=False, name_in_grammar=None, delimiter=""):
        out = False
        for value in self._values.values():
            definition = self._definition
            if definition._save_to_file(
                file, value, always, name_in_grammar, delimiter=delimiter
            ):
                name_in_grammar = False
                delimiter = definition.repeated_delimiter
                out = True
        return out

    def _validate(self, why="save"):
        if why == "save" and not self._definition.is_optional and not self.has_any_value():
            DataValidityError.warn(
                f"Non-optional section {self._definition.name} has no value to save"
            )
        for value in self._values.values():
            value._validate(why)
