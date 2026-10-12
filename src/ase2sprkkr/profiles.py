"""Execution-profile views used by the ase2sprkkr configuration.

Persistent named profiles live in the ase2sprkkr configuration. ``ase`` is a
virtual read-only profile backed by ASE's own configuration. ``legacy`` is a
persisted compatibility profile, but is exposed through a read-only view and
can only be modified through the deprecated configuration options.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
import inspect
import os
import shlex
import warnings

from .common.configuration_transaction import ConfigurationTransaction
from .common.repeated_configuration_containers import RepeatedConfigurationContainer


_PACKAGE_DIR = os.path.normcase(os.path.realpath(os.path.dirname(__file__)))


def _external_warning_stacklevel():
    """Return a warning stacklevel for the first caller outside this package."""
    frame = inspect.currentframe()
    if frame is None:
        return 1
    frame = frame.f_back
    stacklevel = 1
    try:
        while frame is not None:
            filename = os.path.normcase(os.path.realpath(frame.f_code.co_filename))
            try:
                internal = os.path.commonpath((filename, _PACKAGE_DIR)) == _PACKAGE_DIR
            except ValueError:
                internal = False
            if not internal:
                return stacklevel
            frame = frame.f_back
            stacklevel += 1
    finally:
        del frame
    return 1


def ase_command_profile_values(command):
    """Map ASE's canonical ``kkrscf`` command to the common profile fields.

    ASE's ``command`` remains a real command line.  For SPR-KKR it must name
    the canonical serial ``kkrscf`` executable directly.  Its directory and
    suffix locate the other executables from the same installation.  Command
    prefixes and wrappers are intentionally rejected because MPI and execution
    policy can depend on the particular SPR-KKR task.
    """
    if not isinstance(command, str) or not command.strip():
        raise ValueError(
            "ASE SPR-KKR command must directly name the kkrscf executable"
        )

    try:
        argv = shlex.split(command)
    except ValueError as exc:
        raise ValueError("ASE SPR-KKR command is not a valid command line") from exc

    if len(argv) != 1:
        raise ValueError(
            "ASE SPR-KKR command must directly name the kkrscf executable; "
            "command prefixes and wrappers are not supported"
        )

    executable_dir, filename = os.path.split(argv[0])
    if not filename.startswith("kkrscf"):
        raise ValueError(
            "ASE SPR-KKR command must directly name the canonical kkrscf executable"
        )

    values = {"executable_suffix": filename[len("kkrscf"): ]}
    if executable_dir:
        values["executable_dir"] = executable_dir
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
            lambda: profiles.legacy_value("mpi", "auto"), "mpi"
        )
        self.mpi_warning = ReadOnlyProfileValue(
            lambda: profiles.legacy_value("mpi_warning", True), "mpi_warning"
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
            "mpi_warning": self.mpi_warning(),
            "executable_dir": self.executable_dir(),
            "executable_suffix": self.executable_suffix(),
        }

    get = as_dict

    def set(self, *args, **kwargs):
        raise TypeError(
            "The 'legacy' profile is read-only; create a named profile or use "
            "the ASE [sprkkr] profile instead"
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
    legacy_fields = (
        "mpi",
        "mpi_warning",
        "executable_dir",
        "executable_suffix",
    )

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

    def _warn_legacy_write(self):
        if getattr(self, "_legacy_write_warning_emitted", False):
            return
        warnings.warn(
            "config.executables.dir is deprecated together with "
            "config.executables.suffix, config.running.mpi and "
            "config.running.mpi_warning.\n"
            "\n"
            "To migrate the existing configuration to config.yaml and convert "
            "these settings to an execution profile, run:\n"
            "  ase2sprkkr config --migrate\n"
            "\n"
            "For new configuration, use an execution profile directly. For example,\n"
            "configure SPR-KKR through ASE:\n"
            "\n"
            "  [sprkkr]\n"
            "  command = /opt/sprkkr/kkrscf_v8\n"
            "  mpi = 8\n"
            "\n"
            "or create a named ASE2SPRKKR profile:\n"
            "\n"
            "  ase2sprkkr config --set profiles.local "
            "'{mpi: 8, executable_dir: /opt/sprkkr, executable_suffix: _v8}'\n"
            "  ase2sprkkr config --set default_profile local",
            FutureWarning,
            stacklevel=_external_warning_stacklevel(),
        )
        self._legacy_write_warning_emitted = True

    def _set_legacy_value(self, name, value):
        """Stage a deprecated option into the internal legacy profile."""
        self._warn_legacy_write()
        profile = self._legacy_profile(create=True)
        option = getattr(profile, name)
        transaction = ConfigurationTransaction.current(profile)
        if transaction is None:
            option.set(value)
        else:
            option.stage(transaction, value)

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

    @classmethod
    def _reserved_in(cls, values):
        if isinstance(values, str):
            names = (values,)
        elif isinstance(values, Mapping):
            names = values.keys()
        else:
            return ()
        return tuple(sorted(cls.reserved_names.intersection(names)))

    def stage(
        self,
        transaction,
        values=None,
        value=None,
        *,
        unknown="find",
        merge=False,
        **kwargs,
    ):
        """Stage public profile changes without exposing reserved profiles."""
        reserved = self._reserved_in(values)
        if reserved:
            names = ", ".join(repr(name) for name in reserved)
            raise ValueError(f"Reserved profiles are read-only: {names}")

        legacy = self._values.get("legacy")
        changed = super().stage(
            transaction,
            values,
            value,
            unknown=unknown,
            merge=merge,
            **kwargs,
        )
        if legacy is not None and not merge:
            self._values["legacy"] = legacy
            legacy._container = self
        return changed

    def _stage_persisted(
        self,
        transaction,
        values=None,
        value=None,
        *,
        unknown="find",
        merge=False,
        **kwargs,
    ):
        """Stage profiles read from persisted ASE2SPRKKR configuration."""
        return super().stage(
            transaction,
            values,
            value,
            unknown=unknown,
            merge=merge,
            **kwargs,
        )

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
