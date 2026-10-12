"""Execution-profile views used by the ase2sprkkr configuration.

Persistent named profiles live in the ase2sprkkr configuration. ``ase`` and
``legacy`` are virtual read-only profiles. The legacy profile is backed by an
internal normal profile which is written only through the deprecated options;
the ASE profile is backed by ASE's own configuration.
"""

from __future__ import annotations

from collections.abc import Iterator

from .common.configuration_transaction import ConfigurationTransaction
from .common.repeated_configuration_containers import RepeatedConfigurationContainer


class ReadOnlyProfileValue:
    """Small option-like read-only value used by virtual profiles."""

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
    """Read-only public view of the deprecated execution configuration."""

    name = "legacy"
    source = "legacy"
    read_only = True

    def __init__(self, profiles):
        self._profiles = profiles
        self.mpi = ReadOnlyProfileValue(
            lambda: profiles.legacy_value("mpi"), "mpi"
        )
        self.executable_dir = ReadOnlyProfileValue(
            lambda: profiles.legacy_value("executable_dir"), "executable_dir"
        )
        self.executable_suffix = ReadOnlyProfileValue(
            lambda: profiles.legacy_value("executable_suffix"),
            "executable_suffix",
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


class ASEProfile:
    """Read-only view of ASE's ``[sprkkr]`` calculator configuration.

    ASE calculator profiles are intentionally not copied into the ase2sprkkr
    configuration. Keeping this object as a live view also means changes made
    to ASE's in-memory configuration are immediately visible here.
    """

    name = "ase"
    source = "ase"
    read_only = True

    def __init__(self, ase_config):
        self._ase_config = ase_config

    @property
    def _section(self):
        return self._ase_config.parser["sprkkr"]

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name not in self._section:
            raise AttributeError(name)
        return ReadOnlyProfileValue(lambda: self._section.get(name), name)

    @property
    def command(self):
        return ReadOnlyProfileValue(lambda: self._section.get("command"), "command")

    def as_dict(self, *args, **kwargs):
        return dict(self._section.items())

    get = as_dict

    def set(self, *args, **kwargs):
        raise TypeError(
            "The 'ase' profile is read-only in this configuration API; "
            "edit the ASE [sprkkr] configuration instead"
        )


class ProfilesConfigurationContainer(RepeatedConfigurationContainer):
    """Dict-like repeated section with virtual external profiles.

    Normal entries are stored by :class:`RepeatedConfigurationContainer` and
    therefore participate in parsing/serialization. ``legacy`` uses a separate
    internal section so deprecated options can use the same profile schema
    without making that section mutable or persistent through the profiles API.
    """

    reserved_names = frozenset(("legacy", "ase"))
    legacy_fields = ("mpi", "executable_dir", "executable_suffix")

    def _legacy_profile(self, create=False):
        """Return the internal mutable legacy profile used by old options."""
        profile = getattr(self, "_legacy", None)
        if profile is None and create:
            profile = self._definition.create_object(self, repeated=False)
            self._legacy = profile
        return profile

    def _legacy_available(self) -> bool:
        profile = self._legacy_profile()
        return profile is not None and any(
            getattr(profile, name)._value is not None
            for name in self.legacy_fields
        )

    def legacy_value(self, name, default=None):
        """Return a value from the internal legacy profile."""
        profile = self._legacy_profile()
        value = None if profile is None else getattr(profile, name)()
        if value is None and default is not None:
            return default() if callable(default) else default
        return value

    def set_legacy_value(self, name, value):
        """Stage a deprecated option into the internal legacy profile."""
        profile = self._legacy_profile(create=True)
        option = getattr(profile, name)
        transaction = ConfigurationTransaction.current()
        if transaction is None:
            option.set(value)
        else:
            option.stage(transaction, value)

    def legacy_option(self, name, create=False):
        """Return the mutable backing option for compatibility plumbing."""
        profile = self._legacy_profile(create=create)
        return None if profile is None else getattr(profile, name)

    @staticmethod
    def _ase_config():
        try:
            from ase.config import cfg
        except (ImportError, AttributeError):
            return None
        return cfg

    def _ase_available(self) -> bool:
        cfg = self._ase_config()
        if cfg is None:
            return False
        try:
            section = cfg.parser["sprkkr"]
        except KeyError:
            return False
        # ``command`` is the common ASE profile requirement. An incomplete
        # section should not silently become the preferred runnable profile.
        return bool(section.get("command"))

    def _virtual(self, name):
        if name == "ase" and self._ase_available():
            return ASEProfile(self._ase_config())
        if name == "legacy" and self._legacy_available():
            return LegacyProfile(self)
        raise KeyError(name)

    def add(self, id=None):
        if id in self.reserved_names:
            raise ValueError(f"{id!r} is a reserved profile name")
        return super().add(id)

    def __getitem__(self, name):
        try:
            return super().__getitem__(name)
        except (IndexError, KeyError, TypeError):
            return self._virtual(name)

    def __contains__(self, name):
        return super().__contains__(name) or (
            name == "ase" and self._ase_available()
        ) or (
            name == "legacy" and self._legacy_available()
        )

    def _virtual_names(self) -> Iterator[str]:
        if self._ase_available():
            yield "ase"
        if self._legacy_available():
            yield "legacy"

    def __iter__(self):
        # ``items`` and ``values`` intentionally remain inherited and contain
        # persistent entries only. The repeated-container implementation uses
        # them while staging, validating and serializing its contents.
        yield from super().__iter__()
        yield from self._virtual_names()

    def __len__(self):
        return super().__len__() + sum(1 for _ in self._virtual_names())

    def available_items(self):
        """Iterate over both persistent and virtual profiles."""
        for name in self:
            yield name, self[name]

    def resolve_name(self, name=None):
        """Resolve a profile name according to configuration precedence."""
        if name is None:
            name = self._container.default_profile()
        if name is not None:
            if name not in self:
                raise KeyError(f"Unknown execution profile {name!r}")
            return name
        if "ase" in self:
            return "ase"
        if "legacy" in self:
            return "legacy"
        return None

    def resolve(self, name=None):
        """Return the resolved profile object, or ``None`` for auto fallback."""
        name = self.resolve_name(name)
        return None if name is None else self[name]
