"""Execution-profile views used by the ase2sprkkr configuration.

Persistent named profiles live in the ase2sprkkr configuration.  ``ase`` and
``legacy`` are virtual read-only profiles backed by ASE and by the deprecated
ase2sprkkr execution options respectively; they are never serialized as named
ase2sprkkr profiles.
"""

from __future__ import annotations

from collections.abc import Iterator

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


class ASEProfile:
    """Read-only view of ASE's ``[sprkkr]`` calculator configuration.

    ASE calculator profiles are intentionally not copied into the ase2sprkkr
    configuration.  Keeping this object as a live view also means changes made
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
    therefore participate in parsing/serialization.  Virtual entries are added
    only by lookup/iteration, so the base implementation remains the sole
    persistence path and cannot accidentally write them to ase2sprkkr config.
    """

    reserved_names = frozenset(("legacy", "ase"))

    def _legacy_available(self) -> bool:
        config = self._container
        # Explicit presence matters, not whether the value differs from its
        # default.  This also deliberately ignores the suffix environment
        # default: an environment default alone is not a legacy config profile.
        return any(
            getattr(option, "_value", None) is not None
            for option in (
                config.running.mpi,
                config.executables.dir,
                config.executables.suffix,
            )
        )

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
        # ``command`` is the common ASE profile requirement.  An incomplete
        # section should not silently become the preferred runnable profile.
        return bool(section.get("command"))

    def _virtual(self, name):
        if name == "ase" and self._ase_available():
            return ASEProfile(self._ase_config())
        if name == "legacy" and self._legacy_available():
            return LegacyProfile(self._container)
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
        # Iteration is the public list of available profile names.  ``items``
        # and ``values`` intentionally remain inherited and persistent-only,
        # because the base repeated-container implementation uses them while
        # staging, validating and serializing its contents.
        yield from super().__iter__()
        yield from self._virtual_names()

    def __len__(self):
        return super().__len__() + sum(1 for _ in self._virtual_names())

    def available_items(self):
        """Iterate over both persistent and virtual profiles."""
        for name in self:
            yield name, self[name]

    def resolve_name(self, name=None):
        """Resolve the profile name according to configuration precedence.

        An explicit ``name`` wins.  Otherwise an explicitly configured
        ``default_profile`` wins, followed by ASE and then by the deprecated
        legacy configuration.  ``None`` means that no configured profile is
        available and the caller may use its automatic fallback.
        """
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
