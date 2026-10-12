"""Machine-specific SPR-KKR execution profiles.

A runtime profile is deliberately separate from the mutable user configuration:
configuration profiles are resolved once into this small immutable value before
a calculation is prepared.
"""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass, replace
from typing import Any

from ..configuration import config
from ..execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalExecutor,
    PsijExecutor,
    RemoteExecutor,
    SshTransport,
)
from ..profiles import ase_command_profile_values


@dataclass(frozen=True)
class SPRKKRProfile:
    """Describe how and where SPR-KKR executables are invoked."""

    executable_dir: Any = False
    executable_suffix: Any = ""
    mpi: Any = False
    executor: str = "local"
    host: Any = None
    local_root: Any = None
    remote_root: Any = None
    ssh_options: Any = ()
    remote_command: Any = ()
    input_sync: str = "changed"
    output_sync: str = "changed"
    psij_backend: Any = None
    psij_url: Any = None

    def with_overrides(self, **overrides):
        """Return this profile with non-``None`` compatibility overrides applied."""
        allowed = {"mpi", "executable_dir", "executable_suffix"}
        unknown = set(overrides) - allowed
        if unknown:
            names = ", ".join(sorted(unknown))
            raise TypeError(f"Unknown SPR-KKR profile override(s): {names}")

        if "mpi" in overrides and overrides["mpi"] is not None:
            overrides["mpi"] = _normalize_mpi(overrides["mpi"])
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
        cache_files=(),
        mpi_supported=False,
        gdb=False,
    ):
        """Create an execution specification for one SPR-KKR task.

        The profile describes MPI intent only. The concrete launcher is chosen
        by the executor on the machine where the calculation actually runs.
        """
        executable += self.executable_suffix or ""
        if self.executable_dir:
            executable = os.path.join(self.executable_dir, executable)

        launch = LaunchSpec.from_mpi(self.mpi) if mpi_supported else LaunchSpec()
        argv = [executable]
        mpi_argv = None
        stdin = input_file

        if launch:
            input_name = input_file.name
            # SPR-KKR limits the input-file command-line argument to 80 chars.
            if len(input_name) > 80:
                input_name = os.path.relpath(input_name, directory)
            mpi_argv = [executable + "MPI", input_name]

        if gdb:
            if launch:
                raise ValueError("gdb with deferred MPI execution is not supported")
            stdin = None
            print(f"run <{input_file.name}")
            argv = ["gdb", executable]

        return ExecutionSpec(
            argv=tuple(argv),
            input_file=input_file,
            output_file=output_file,
            directory=directory,
            stdin=stdin,
            input_files=tuple(input_files),
            output_files=tuple(output_files),
            cache_files=tuple(cache_files),
            mpi_argv=None if mpi_argv is None else tuple(mpi_argv),
            launch=launch,
        )


def _profile_value(profile, name):
    """Read one option-like value from a configuration profile."""
    try:
        value = getattr(profile, name)
    except AttributeError:
        return None
    return value() if callable(value) else value


def _shell_argv(value):
    """Normalize a shell-style command string or argv sequence to a tuple."""
    if not value:
        return ()
    if isinstance(value, str):
        return tuple(shlex.split(value))
    return tuple(value)


def _normalize_mpi(value):
    """Normalize an explicit MPI command while retaining symbolic modes."""
    if isinstance(value, str) and value != "auto":
        return _shell_argv(value)
    return value


def profile_from_ase_command(command):
    """Translate ASE's ``[sprkkr] command`` into an SPR-KKR profile."""
    return SPRKKRProfile(**ase_command_profile_values(command))


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
        mpi = False if profile.mpi is None else _normalize_mpi(profile.mpi)
        if mpi != profile.mpi:
            changes["mpi"] = mpi
        if not profile.executor:
            changes["executor"] = "local"
        ssh_options = _shell_argv(profile.ssh_options)
        if ssh_options != profile.ssh_options:
            changes["ssh_options"] = ssh_options
        remote_command = _shell_argv(profile.remote_command)
        if remote_command != profile.remote_command:
            changes["remote_command"] = remote_command
        if not profile.input_sync:
            changes["input_sync"] = "changed"
        if not profile.output_sync:
            changes["output_sync"] = "changed"
        return replace(profile, **changes) if changes else profile

    if getattr(profile, "source", None) == "ase":
        profile = profile.configuration_profile()

    suffix = _profile_value(profile, "executable_suffix")
    if suffix is None:
        suffix = os.environ.get("SPRKKR_EXECUTABLE_SUFFIX", "")
    executable_dir = _profile_value(profile, "executable_dir")
    if executable_dir is None:
        executable_dir = False
    mpi = _profile_value(profile, "mpi")
    if mpi is None:
        mpi = False
    else:
        mpi = _normalize_mpi(mpi)

    return SPRKKRProfile(
        executable_dir=executable_dir,
        executable_suffix=suffix,
        mpi=mpi,
        executor=_profile_value(profile, "executor") or "local",
        host=_profile_value(profile, "host"),
        local_root=_profile_value(profile, "local_root"),
        remote_root=_profile_value(profile, "remote_root"),
        ssh_options=_shell_argv(_profile_value(profile, "ssh_options")),
        remote_command=_shell_argv(_profile_value(profile, "remote_command")),
        input_sync=_profile_value(profile, "input_sync") or "changed",
        output_sync=_profile_value(profile, "output_sync") or "changed",
        psij_backend=_profile_value(profile, "psij_backend"),
        psij_url=_profile_value(profile, "psij_url"),
    )


def executor_from_profile(profile: SPRKKRProfile):
    """Create an execution backend from a resolved runtime profile."""
    if profile.executor == "local":
        return LocalExecutor()

    if profile.executor == "ssh":
        if not profile.host:
            raise ValueError("SSH execution profile requires 'host'")
        kwargs = {
            "local_root": profile.local_root,
            "remote_root": profile.remote_root,
            "input_files": profile.input_sync,
            "output_files": profile.output_sync,
        }
        if profile.remote_command:
            kwargs["command"] = profile.remote_command
        if profile.psij_backend is not None or profile.psij_url is not None:
            kwargs["backend"] = "psij"
            kwargs["backend_options"] = {
                "backend": profile.psij_backend or "slurm",
                "url": profile.psij_url,
            }
        return RemoteExecutor(
            SshTransport(profile.host, options=profile.ssh_options),
            **kwargs,
        )

    if profile.executor == "psij":
        return PsijExecutor(profile.psij_backend or "slurm", url=profile.psij_url)

    raise ValueError(f"Unknown execution backend {profile.executor!r}")


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
