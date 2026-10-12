"""Machine-specific SPR-KKR execution profiles.

A runtime profile is deliberately separate from the mutable user configuration:
configuration profiles are resolved once into this small immutable value before
a calculation is prepared.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, replace
from typing import Any

from ..configuration import config
from ..execution import ExecutionSpec


@dataclass(frozen=True)
class SPRKKRProfile:
    """Describe how SPR-KKR executables are invoked on one computer."""

    executable_dir: Any = False
    executable_suffix: Any = ""
    mpi: Any = False

    def with_overrides(self, **overrides):
        """Return this profile with non-``None`` overrides applied.

        The profile is immutable, so actual changes create a copy. If all
        supplied values are ``None`` or already equal to the current values,
        return ``self`` instead. ``None`` therefore means "keep the profile
        value", which matches the compatibility arguments of ``SPRKKR``.
        """
        allowed = {"mpi", "executable_dir", "executable_suffix"}
        unknown = set(overrides) - allowed
        if unknown:
            names = ", ".join(sorted(unknown))
            raise TypeError(f"Unknown SPR-KKR profile override(s): {names}")

        changes = {
            name: value
            for name, value in overrides.items()
            if value is not None and value != getattr(self, name)
        }
        return replace(self, **changes) if changes else self

    def create_execution(
        self,
        executable,
        input_file,
        output_file,
        directory,
        *,
        input_files=(),
        output_files=(),
        mpi_supported=False,
        gdb=False,
    ):
        """Create a concrete execution specification for one SPR-KKR task.

        ``mpi`` stored in this profile must already be resolved to a concrete
        command by the executor. Resolving ``True``, ``"auto"`` or a process
        count belongs to the executor, because only the executor knows the
        machine on which the command will run.
        """
        executable += self.executable_suffix or ""
        if self.executable_dir:
            executable = os.path.join(self.executable_dir, executable)

        mpi_command = self.mpi if mpi_supported else None

        if mpi_command:
            input_name = input_file.name
            # SPR-KKR limits the input-file command-line argument to 80 chars.
            if len(input_name) > 80:
                input_name = os.path.relpath(input_name, directory)
            argv = [*mpi_command, executable + "MPI", input_name]
            stdin = None
        else:
            argv = [executable]
            stdin = input_file

        if gdb:
            stdin = None
            if mpi_command:
                print(f"run {input_file.name}")
                argv = argv[:-1]
            else:
                print(f"run <{input_file.name}")
            argv = ["gdb", *argv]

        return ExecutionSpec(
            argv=tuple(argv),
            input_file=input_file,
            output_file=output_file,
            directory=directory,
            stdin=stdin,
            input_files=tuple(input_files),
            output_files=tuple(output_files),
        )


def _profile_value(profile, name):
    """Read one option-like value from a configuration profile."""
    try:
        value = getattr(profile, name)
    except AttributeError:
        return None
    return value() if callable(value) else value


def profile_from_ase_command(command):
    """Translate ASE's ``[sprkkr] command`` into an SPR-KKR profile.

    A command without ``{}`` denotes the directory containing the SPR-KKR
    executables. A command of the form ``directory/{}suffix`` additionally
    specifies the executable suffix. The placeholder is intentionally much
    more restricted than a general format string: it must occur exactly once
    and at the beginning of the final path component.
    """
    if not isinstance(command, str) or not command:
        raise ValueError("ASE SPR-KKR command must be a non-empty path")

    count = command.count("{}")
    if count == 0:
        return SPRKKRProfile(executable_dir=command)
    if count != 1:
        raise ValueError("ASE SPR-KKR command may contain at most one '{}' placeholder")

    executable_dir, filename = os.path.split(command)
    if "{}" in executable_dir or not filename.startswith("{}"):
        raise ValueError(
            "ASE SPR-KKR command placeholder must start the final path component, "
            "for example '/opt/sprkkr/{}.mpi'"
        )
    return SPRKKRProfile(
        executable_dir=executable_dir or False,
        executable_suffix=filename[2:],
    )


def _configured_profile(profile, configuration):
    if profile is None or isinstance(profile, str):
        profile = configuration.profiles.resolve(profile)

    if profile is None:
        return SPRKKRProfile(
            executable_suffix=os.environ.get("SPRKKR_EXECUTABLE_SUFFIX", ""),
        )

    if isinstance(profile, SPRKKRProfile):
        changes = {}
        if profile.executable_dir is None:
            changes["executable_dir"] = False
        if profile.executable_suffix is None:
            changes["executable_suffix"] = ""
        if profile.mpi is None:
            changes["mpi"] = False
        return replace(profile, **changes) if changes else profile

    if getattr(profile, "source", None) == "ase":
        return profile_from_ase_command(profile.command())

    suffix = _profile_value(profile, "executable_suffix")
    if suffix is None:
        suffix = os.environ.get("SPRKKR_EXECUTABLE_SUFFIX", "")
    executable_dir = _profile_value(profile, "executable_dir")
    if executable_dir is None:
        executable_dir = False
    mpi = _profile_value(profile, "mpi")
    if mpi is None:
        mpi = False

    return SPRKKRProfile(
        executable_dir=executable_dir,
        executable_suffix=suffix,
        mpi=mpi,
    )


def resolve_sprkkr_profile(profile=None, *, configuration=None):
    """Resolve a user/configuration profile into an immutable runtime profile.

    Resolution order for ``profile=None`` is the one defined by the profile
    configuration container: configured ``default_profile``, ASE's virtual
    profile, the deprecated legacy profile, then built-in defaults. A string
    explicitly selects a named profile. Per-calculation compatibility
    overrides are applied by ``SPRKKR`` after this resolution.
    """
    configuration = configuration or config
    return _configured_profile(profile, configuration)
