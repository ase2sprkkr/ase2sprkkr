"""Execution-profile views used by the ase2sprkkr configuration.

Persistent named profiles live in the ase2sprkkr configuration. ``ase`` is a
virtual read-only profile backed by ASE's own configuration. ``legacy`` is a
persisted compatibility profile, but is exposed through a read-only view and
can only be modified through the deprecated configuration options.
"""

from __future__ import annotations

from collections.abc import Iterator
import os

from .common.configuration_transaction import ConfigurationTransaction
from .common.repeated_configuration_containers import RepeatedConfigurationContainer


def ase_command_profile_values(command):
    """Map ASE's ``[sprkkr] command`` to the common profile fields.

    A command without ``{}`` denotes the executable directory. A command of
    the form ``directory/{}suffix`` additionally specifies the executable
    suffix. The placeholder is intentionally restricted to the beginning of
    the final path component; it is not a general formatting language.
    """
    if not isinstance(command, str) or not command:
        raise ValueError("ASE SPR-KKR command must be a non-empty path")

    count = command.count("{}")
    if count == 0:
        return {"executable_dir": command}
    if count != 1:
        raise ValueError("ASE SPR-KKR command may contain at most one '{}' placeholder")

    executable_dir, filename = os.path.split(command)
    if "{}" in executable_dir or not filename.startswith("{}"):
        raise ValueError(
            "ASE SPR-KKR command placeholder must start the final path component, "
            "for example '/opt/sprkkr/{}.mpi'"
        )

    values = {}
    if executable_dir:
        values["executable_dir"] = executable_dir
    suffix = filename[2:]
    if suffix:
        values["executable_suffix"] = suffix
    return values


class ReadOnlyProfileValue:
    """Small option-like read-only value used by virtual/read-only profiles."""

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
    """Typed read-only view of ASE's ``[sprkkr]`` calculator configuration."""

    name = "ase"
    source = "ase"
    read_only = True

    def __init__(self, profiles, ase_config):
        self._profiles = profiles
        self._ase_config = ase_config

    @property
    def _section(self):
        return self._ase_config.parser["sprkkr"]

    def configuration_profile(self):
        """Return a typed, non-persistent profile snapshot."""
        return self._profiles._parse_ase_profile(self._ase_config)

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            getattr(self.configuration_profile(), name)
        except AttributeError:
            raise AttributeError(name) from None
        return ReadOnlyProfileValue(
            lambda: getattr(self.configuration_profile(), name)(),
            name,
        )

    @property
    def command(self):
        return ReadOnlyProfileValue(lambda: self._section.get("command"), "command")

    def as_dict(self, *args, **kwargs):
        values = self.configuration_profile().as_dict(
            only_changed="explicit", generated=False
        ) or {}
        return {"command": self.command(), **values}

    get = as_dict

    def set(self, *args, **kwargs):
        raise TypeError(
            "The 'ase' profile is read-only in this configuration API; "
            "edit the ASE [sprkkr] configuration instead"
        )


class ProfilesConfigurationContainer(RepeatedConfigurationContainer):
    """Dict-like repeated section with persistent and virtual profile views."""

    reserved_names = frozenset(("legacy", "ase"))
    legacy_fields = ("mpi", "executable_dir", "executable_suffix")

    def _add(self, id=None):
        """Create one stored profile while keeping ``ase`` strictly virtual."""
        if id == "ase":
            raise ValueError("'ase' is a virtual profile and cannot be stored")
        return super()._add(id)

    def _legacy_profile(self, create=False):
        """Return the mutable backing profile used by deprecated options."""
        profile = self._values.get("legacy")
        if profile is None and create:
            profile = self._add("legacy")
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
        transaction = ConfigurationTransaction.current(profile)
        if transaction is None:
            option.set(value)
        else:
            option.stage(transaction, value)

    def legacy_option(self, name, create=False):
        """Return mutable backing option for compatibility plumbing only."""
        profile = self._legacy_profile(create=create)
        return None if profile is None else getattr(profile, name)

    @staticmethod
    def _ase_config():
        try:
            from ase.config import cfg
        except (ImportError, AttributeError):
            return None
        return cfg

    def _parse_ase_profile(self, cfg):
        """Parse ASE's text configuration through the common profile schema."""
        section = cfg.parser["sprkkr"]
        command = section.get("command")
        values = dict(section.items())
        values.pop("command", None)

        duplicate = {"executable_dir", "executable_suffix"} & values.keys()
        if duplicate:
            names = ", ".join(sorted(duplicate))
            raise ValueError(
                f"ASE [sprkkr] {names} must be expressed through the command option"
            )

        profile = self._definition.create_object(self, repeated=False)
        parsed = {}
        for name, value in values.items():
            try:
                option = getattr(profile, name)
            except AttributeError:
                raise KeyError(
                    f"Unknown ASE SPR-KKR profile option {name!r}"
                ) from None
            parsed[name] = option._definition.type.parse(value)

        parsed.update(ase_command_profile_values(command))
        profile.set(parsed, unknown="fail")
        return profile

    def _ase_available(self) -> bool:
        cfg = self._ase_config()
        if cfg is None:
            return False
        try:
            section = cfg.parser["sprkkr"]
        except KeyError:
            return False
        return bool(section.get("command"))

    def add(self, id=None):
        if id in self.reserved_names:
            raise ValueError(f"{id!r} is a reserved profile name")
        return self._add(id)

    def __getitem__(self, name):
        if name == "legacy":
            if self._legacy_available():
                return LegacyProfile(self)
            raise KeyError(name)
        if name == "ase":
            if self._ase_available():
                return ASEProfile(self, self._ase_config())
            raise KeyError(name)
        return super().__getitem__(name)

    def __contains__(self, name):
        if name == "legacy":
            return self._legacy_available()
        if name == "ase":
            return self._ase_available()
        return super().__contains__(name)

    def _virtual_names(self) -> Iterator[str]:
        if self._ase_available():
            yield "ase"

    def __iter__(self):
        for name in super().__iter__():
            if name != "legacy" or self._legacy_available():
                yield name
        yield from self._virtual_names()

    def __len__(self):
        return sum(1 for _ in self)

    def items(self):
        """Iterate over public profile views only."""
        for name in self:
            yield name, self[name]

    def values(self):
        """Iterate over public profile views only."""
        for _name, profile in self.items():
            yield profile

    def available_items(self):
        """Backward-compatible alias for public profile items."""
        yield from self.items()

    def _get_members(
        self,
        name,
        unknown,
        is_option,
        lower_case,
        *,
        accept,
    ):
        # The backing legacy section is an implementation detail. Generic
        # lookup must not expose mutable options below its read-only public view.
        head = name.split(".", 1)[0] if isinstance(name, str) else name
        if head == "legacy":
            return
        yield from super()._get_members(
            name,
            unknown,
            is_option,
            lower_case,
            accept=accept,
        )

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
