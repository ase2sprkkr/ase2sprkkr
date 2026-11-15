"""Execution-profile views used by the ase2sprkkr configuration.

Persistent named profiles live in the ase2sprkkr configuration. ``default``
is a virtual alias for the selected execution profile. ``ase`` is a
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

import pyparsing as pp

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


class ASEProfileValue(ReadOnlyProfileValue):
    """Option-like field backed by ASE's own INI configuration."""

    def __init__(self, profile, name, getter):
        super().__init__(getter, name)
        self._profile = profile

    def set(self, value):
        """Validate and stage a value without writing to disk."""
        self._profile.set({self.name: value})

    def clear(self):
        """Stage removal of a value without writing to disk."""
        self._profile.unset(self.name)

    def store_permanent(self, *, system=False):
        """Save pending ASE configuration changes, permitting /etc if requested."""
        self._profile.store_permanent(system=system)

    def set_permanent(self, value, *, system=False):
        """Set and save a value, optionally permitting system-file writes."""
        self.set(value)
        self.store_permanent(system=system)


class ASEProfile:
    """Editable, virtual view of ASE's [sprkkr] configuration.

    Values are edited in memory and saved explicitly to ASE's own
    configuration file using store_permanent() or the root
    save_user_preferences(). The ASE profile is never stored in YAML.
    """

    name = "ase"
    source = "ase"
    read_only = False

    def __init__(self, profiles, ase_config):
        self._profiles = profiles
        self._ase_config = ase_config

    @property
    def _section(self):
        parser = self._ase_config.parser
        return parser["sprkkr"] if parser.has_section("sprkkr") else None

    def configuration_profile(self):
        """Return a typed view without changing the ASE configuration."""
        if self._profiles._ase_available():
            return self._profiles._parse_ase_profile(self._ase_config)
        return self._profiles._definition.create_object(
            self._profiles, repeated=False
        )

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            getattr(self.configuration_profile(), name)
        except AttributeError:
            raise AttributeError(name) from None
        return ASEProfileValue(
            self, name, lambda: getattr(self.configuration_profile(), name)(),
        )

    @property
    def command(self):
        return ASEProfileValue(
            self, "command",
            lambda: self._section.get("command") if self._section else None,
        )

    def as_dict(self, *args, **kwargs):
        values = self.configuration_profile().as_dict(
            only_changed="explicit", generated=False
        ) or {}
        return {"command": self.command(), **values}

    get = as_dict

    @staticmethod
    def _ini_value(value):
        """Represent a typed profile value for ASE's INI parser."""
        if isinstance(value, bool):
            return "true" if value else "false"
        return str(value)

    def _edit(self, changes=None, *, remove=()):
        """Validate a proposed ASE profile and stage the change in memory."""
        import copy
        import shlex
        from types import SimpleNamespace

        changes = dict(changes or {})
        removals = tuple(remove)
        if not changes and not removals:
            return
        if changes.keys() & set(removals):
            raise ValueError("Cannot set and unset the same ASE profile option")

        for name in (*changes, *removals):
            if name == "command":
                continue
            try:
                getattr(self.configuration_profile(), name)
            except AttributeError:
                raise ValueError(
                    f"Unknown ASE execution-profile option: {name}"
                ) from None

        if "command" in changes and (
            {"executable_dir", "executable_suffix"} & changes.keys()
        ):
            raise ValueError(
                "Specify either command or executable_dir/executable_suffix"
            )
        if "command" in removals:
            raise ValueError(
                "Unsetting the ASE command would leave the profile unusable; "
                "use another command or edit ASE config.ini directly"
            )

        original_command = self.command()
        original_parts = (
            ase_command_profile_values(original_command)
            if original_command else {"executable_dir": "", "executable_suffix": ""}
        )
        updates = {}
        if "command" in changes:
            command = changes.pop("command")
            ase_command_profile_values(command)
            updates["command"] = command
        elif (
            {"executable_dir", "executable_suffix"} &
            (changes.keys() | set(removals))
        ) or not original_command:
            directory = changes.pop(
                "executable_dir", original_parts.get("executable_dir", "")
            )
            suffix = changes.pop(
                "executable_suffix", original_parts.get("executable_suffix", "")
            )
            if "executable_dir" in removals:
                directory = ""
            if "executable_suffix" in removals:
                suffix = ""
            if not isinstance(directory, str) or not isinstance(suffix, str):
                raise TypeError("Executable directory and suffix must be strings")
            filename = "kkrscf" + suffix
            filename = os.path.join(directory, filename) if directory else filename
            command = shlex.quote(filename)
            ase_command_profile_values(command)
            updates["command"] = command

        for name, value in changes.items():
            updates[name] = self._ini_value(value)
        for name in removals:
            if name not in ("executable_dir", "executable_suffix"):
                updates[name] = None

        for value in updates.values():
            if value is not None and any(char in value for char in (chr(10), chr(13))):
                raise ValueError("ASE profile INI values must be single-line")

        parser = copy.deepcopy(self._ase_config.parser)
        if not parser.has_section("sprkkr"):
            parser.add_section("sprkkr")
        section = parser["sprkkr"]
        for name, value in updates.items():
            if value is None:
                section.pop(name, None)
            else:
                section[name] = value

        # One shared schema validates command syntax and profile options.
        self._profiles._parse_ase_profile(SimpleNamespace(parser=parser))

        # Stage only after successful validation: rejected changes do not
        # affect the parser or pending writes. All ASEProfile views share the
        # same staged state through the profiles container.
        pending = getattr(self._profiles, "_ase_pending_updates", {})
        self._profiles._ase_pending_updates = {**pending, **updates}
        self._ase_config.parser = parser

    def set(self, values=None, **kwargs):
        """Validate and stage ASE profile fields without saving to disk."""
        from collections.abc import Mapping
        if values is None:
            values = {}
        if not isinstance(values, Mapping):
            raise TypeError("ASE profile values must be a mapping")
        proposed = dict(values)
        proposed.update(kwargs)
        self._edit(proposed)

    def unset(self, name):
        """Stage removal of one ASE profile field without saving."""
        self._edit(remove=(name,))

    def store_permanent(self, *, system=False):
        """Persist all staged ASE fields to ASE's authoritative INI file."""
        self._profiles.save_ase_changes(system=system)

    def set_permanent(self, values=None, *, system=False, **kwargs):
        """Set and persist ASE profile fields in one operation."""
        self.set(values, **kwargs)
        self.store_permanent(system=system)


class ProfilesConfigurationContainer(RepeatedConfigurationContainer):
    """Dict-like repeated section with persistent and virtual profile views."""

    class RemovalChange:
        """Rollback state for removing one stored profile."""

        def __init__(self, profiles):
            self.profiles = profiles
            self.old_values = profiles._values.copy()

        def commit(self):
            return None

        def rollback(self):
            self.profiles._values = self.old_values

    reserved_names = frozenset(("legacy", "ase", "default"))
    legacy_fields = (
        "mpi",
        "mpi_warning",
        "executable_dir",
        "executable_suffix",
    )

    def _add(self, id=None):
        """Create one stored profile while keeping ``ase`` strictly virtual."""
        if id in ("ase", "default"):
            raise ValueError(f"{id!r} is a virtual profile and cannot be stored")
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

    @staticmethod
    def _parse_ase_value(section, name, option, value):
        """Parse one ASE INI value using the common profile type contract."""
        value_type = option._definition.type
        try:
            return value_type.parse(value)
        except pp.ParseBaseException as exc:
            try:
                return section.getboolean(name)
            except ValueError:
                raise ValueError(
                    f"Invalid ASE SPR-KKR profile option {name!r}: {value!r}"
                ) from exc

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
            parsed[name] = self._parse_ase_value(section, name, option, value)

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
            raise ValueError(f"Reserved profiles are read-only in ASE2SPRKKR YAML: {names}. Edit ASE settings via profiles.ase.")

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

    def remove(self, name):
        """Transactionally remove one mutable named profile."""
        if name in self.reserved_names:
            raise ValueError(f"{name!r} is a reserved profile name")
        if name not in self._values:
            raise KeyError(name)
        if self._container.default_profile() == name:
            raise ValueError(
                f"Cannot remove default execution profile {name!r}; "
                "unset or change default_profile first"
            )

        with self._mutation("set") as (transaction, _policy):
            change = self.RemovalChange(self)
            del self._values[name]
            transaction.push(change)

    @property
    def ase(self):
        """ASE-owned profile, editable via the ASE configuration file."""
        return self["ase"]

    @property
    def default(self):
        """Selected profile; a virtual alias, not a stored configuration."""
        return self["default"]

    def __getitem__(self, name):
        if name == "default":
            # This is a virtual name, not a serialized profile.
            return self[self.resolve_name() or "ase"]
        if name == "legacy":
            if self._legacy_available():
                return LegacyProfile(self)
            raise KeyError(name)
        if name == "ase":
            cfg = self._ase_config()
            if cfg is not None:
                return ASEProfile(self, cfg)
            raise KeyError("ASE configuration support is not installed")
        return super().__getitem__(name)

    def __contains__(self, name):
        if name == "default":
            return self.resolve_name() is not None or self._ase_config() is not None
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

    def _as_dict(self, get):
        """Serialize stored profiles while preserving empty named profiles."""
        out = {}
        for name, profile in self._raw_items():
            value = profile._as_dict(get)
            if value is not None:
                out[name] = value
            elif name not in self.reserved_names:
                out[name] = {}
        return out or None

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
        if head in ("legacy", "default"):
            return
        yield from super()._get_members(
            name,
            unknown,
            is_option,
            lower_case,
            accept=accept,
        )

    def resolve_name(self, name=None):
        """Resolve a profile name according to configuration precedence.

        The reserved name "default" always means the currently selected
        profile, never a persistent profile called "default".
        """
        if name in (None, "default"):
            name = self._container.default_profile()
            if name == "default":
                raise ValueError("default_profile cannot refer to the 'default' alias")
        if name is not None:
            if name not in self:
                raise KeyError(f"Unknown execution profile {name!r}")
            return name
        if "ase" in self:
            return "ase"
        if "legacy" in self:
            return "legacy"
        return None

    def default_write_target(self):
        """Resolve a write to the default profile without redundant state.

        An explicitly selected named profile remains the write target.
        Otherwise ASE's [sprkkr] settings are the natural user default.
        A write to a missing ASE profile creates that section in the
        authoritative ASE config.ini, rather than a parallel YAML copy.
        """
        name = self._container.default_profile()
        if name == "default":
            raise ValueError("default_profile cannot refer to the default alias")
        if name == "legacy":
            raise ValueError(
                "The legacy profile is read-only; select a named profile "
                "or edit profiles.ase instead"
            )
        return name or "ase"

    def save_ase_changes(self, *, system=False):
        """Write staged ASE changes, leaving the staging intact on failure.

        This is the single persistence boundary for ASE-owned profiles,
        shared by explicit profile saving, CLI writes, and the root
        Configuration.save_user_preferences() method.
        """
        from .ase_profile_storage import write_ase_sprkkr

        updates = getattr(self, "_ase_pending_updates", {})
        if not updates:
            return
        write_ase_sprkkr(updates, system=system)
        self._ase_pending_updates = {}

    def resolve(self, name=None):
        """Return the resolved profile object, or ``None`` for auto fallback."""
        name = self.resolve_name(name)
        return None if name is None else self[name]
