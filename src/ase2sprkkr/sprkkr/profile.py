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
from kkrunner import ExecutionSpec
from kkrunner.backends.config import BackendConfig
from kkrunner.client import RunClient
from kkrunner.client.pool import default_client_pool
from kkrunner.client.staging import StagingSpec
from kkrunner.client.transports import LocalTransport, SshTransport
from kkrunner.execution import MpiSpec
from ..profiles import ase_command_profile_values


@dataclass(frozen=True)
class LaunchParameters:
    """Profile-level parameters controlling how calculations are launched."""

    executable_dir: Any = False
    executable_suffix: Any = ""
    mpi: Any = "auto"
    local_root: Any = None
    remote_root: Any = None
    input_sync: str = "changed"
    output_sync: str = "changed"
    mpi_warning: bool = True


@dataclass(frozen=True)
class SPRKKRProfile:
    """Describe how and where SPR-KKR executables are invoked.

    The public/runtime profile stays flat for configuration ergonomics. Code
    that creates execution objects splits it into :attr:`transport`, :attr:`backend_config`, and
    per-run staging and MPI options internally.
    """

    executable_dir: Any = False
    executable_suffix: Any = ""
    mpi: Any = "auto"
    executor: str = "local"
    backend: str | None = None
    host: Any = None
    local_root: Any = None
    remote_root: Any = None
    ssh_options: Any = ()
    remote_command: Any = ()
    remote_python: Any = ("python3",)
    input_sync: str = "changed"
    output_sync: str = "changed"
    psij_backend: Any = None
    psij_url: Any = None
    mpi_warning: bool = True
    ssh_zstd: bool = False
    filesystem_id: str | None = None

    @property
    def transport(self):
        """Client connection, independent of the selected execution backend."""
        if self.executor == "ssh" or (self.executor == "psij" and self.host):
            if not self.host:
                raise ValueError("SSH execution profile requires 'host'")
            python_command = self.remote_python
            # The old kkrunner remote CLI was replaced by Python RPC in wip2.
            if self.remote_command and self.remote_command != ("kkrunner", "remote"):
                raise ValueError(
                    "remote_command is no longer supported by KKRunner; "
                    "configure remote_python with the remote Python interpreter instead"
                )
            return SshTransport(
                host=self.host,
                options=self.ssh_options,
                zstd=self.ssh_zstd,
                python_command=python_command,
            )
        if self.executor in ("local", "psij"):
            return LocalTransport()
        raise ValueError(f"Unknown execution location {self.executor!r}")

    @property
    def backend_config(self):
        """Backend identity is separate from the client transport and staging."""
        name = self.backend or (
            "psij" if self.executor == "psij" else
            "process" if self.executor == "ssh" else "direct"
        )
        if name not in ("direct", "process", "psij"):
            raise ValueError(f"Unknown KKRunner backend {name!r}")
        options = (
            {"executor_name": self.psij_backend or "slurm", "url": self.psij_url}
            if name == "psij" else {}
        )
        return BackendConfig(name, options)

    @property
    def launch_parameters(self):
        return LaunchParameters(
            executable_dir=self.executable_dir,
            executable_suffix=self.executable_suffix,
            mpi=self.mpi,
            local_root=self.local_root,
            remote_root=self.remote_root,
            input_sync=self.input_sync,
            output_sync=self.output_sync,
            mpi_warning=self.mpi_warning,
        )

    def with_overrides(self, **overrides):
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
        error_file=None,
        input_files=(),
        output_files=(),
        cache_files=(),
        mpi_supported=False,
        gdb=False,
    ):
        """Create the concrete execution specification for one run."""
        executable += self.executable_suffix or ""
        if self.executable_dir:
            executable = os.path.join(self.executable_dir, executable)

        mpi = MpiSpec.from_value(self.mpi) if mpi_supported else MpiSpec()
        if mpi_supported:
            mpi = replace(mpi, warning=self.mpi_warning)
        staging = StagingSpec(
            local_root=self.local_root,
            remote_root=self.remote_root,
            input_sync=self.input_sync,
            output_sync=self.output_sync,
            filesystem_id=self.filesystem_id,
        )
        argv = [executable]
        mpi_argv = None
        stdin = input_file

        if mpi.requested:
            input_name = input_file.name
            if len(input_name) > 80:
                input_name = os.path.relpath(input_name, directory)
            mpi_argv = [executable + "MPI", input_name]
            if gdb:
                # Let the executor prepend whichever MPI launcher it uses.
                # Each rank gets its own interactive GDB terminal.
                mpi_argv = ["xterm", "-e", "gdb", "--args", *mpi_argv]

        if gdb:
            # This is also the serial command used by ``mpi='auto'`` when MPI
            # autodetection falls back to non-MPI execution.
            stdin = None
            print(f"run <{input_file.name}")
            argv = ["gdb", executable]

        return ExecutionSpec(
            argv=tuple(argv),
            input_file=input_file,
            output_file=output_file,
            error_file=error_file,
            directory=directory,
            stdin=stdin,
            input_files=tuple(input_files),
            output_files=tuple(output_files),
            cache_files=tuple(cache_files),
            mpi_argv=None if mpi_argv is None else tuple(mpi_argv),
            mpi=mpi,
            staging=staging,
            mpi_uses_stdin=False if mpi_argv is not None else True,
        )


def _profile_value(profile, name):
    try:
        value = getattr(profile, name)
    except AttributeError:
        return None
    return value() if callable(value) else value


def _shell_argv(value):
    if value is None:
        return ()
    if isinstance(value, str):
        return tuple(shlex.split(value))
    return tuple(value)


def _normalize_mpi(value):
    if isinstance(value, MpiSpec):
        return value
    if isinstance(value, str):
        return value if value == "auto" else _shell_argv(value)
    if value is None or isinstance(value, (bool, int)):
        return value
    try:
        return tuple(value)
    except TypeError:
        return value


def profile_from_ase_command(command):
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
        mpi = "auto" if profile.mpi is None else _normalize_mpi(profile.mpi)
        if type(mpi) is not type(profile.mpi) or mpi != profile.mpi:
            changes["mpi"] = mpi
        if profile.mpi_warning is None:
            changes["mpi_warning"] = True
        if not profile.executor:
            changes["executor"] = "local"
        remote_python = _shell_argv(profile.remote_python) or ("python3",)
        if remote_python != profile.remote_python:
            changes["remote_python"] = remote_python
        ssh_options = _shell_argv(profile.ssh_options)
        if type(ssh_options) is not type(profile.ssh_options) or ssh_options != profile.ssh_options:
            changes["ssh_options"] = ssh_options
        remote_command = _shell_argv(profile.remote_command)
        if type(remote_command) is not type(profile.remote_command) or remote_command != profile.remote_command:
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
        mpi = "auto"
    else:
        mpi = _normalize_mpi(mpi)
    mpi_warning = _profile_value(profile, "mpi_warning")
    if mpi_warning is None:
        mpi_warning = True

    return SPRKKRProfile(
        executable_dir=executable_dir,
        executable_suffix=suffix,
        mpi=mpi,
        executor=_profile_value(profile, "executor") or "local",
        backend=_profile_value(profile, "backend"),
        host=_profile_value(profile, "host"),
        local_root=_profile_value(profile, "local_root"),
        remote_root=_profile_value(profile, "remote_root"),
        ssh_options=_shell_argv(_profile_value(profile, "ssh_options")),
        ssh_zstd=bool(_profile_value(profile, "ssh_zstd")),
        remote_command=_shell_argv(_profile_value(profile, "remote_command")),
        remote_python=_shell_argv(_profile_value(profile, "remote_python")) or ("python3",),
        input_sync=_profile_value(profile, "input_sync") or "changed",
        output_sync=_profile_value(profile, "output_sync") or "changed",
        psij_backend=_profile_value(profile, "psij_backend"),
        psij_url=_profile_value(profile, "psij_url"),
        filesystem_id=_profile_value(profile, "filesystem_id"),
        mpi_warning=mpi_warning,
    )


def create_client_from_profile(profile: SPRKKRProfile):
    """Create a client combining independent transport and backend configuration."""
    return RunClient.from_config(profile.transport, profile.backend_config)


def client_from_profile(profile: SPRKKRProfile):
    """Reuse a client based on its backend and connection, not staging or MPI."""
    return default_client_pool().for_config(profile.transport, profile.backend_config)


# Retain the calculator's existing executor= argument and import paths while
# treating its value as a KKRunner RunClient from now on.
create_executor_from_profile = create_client_from_profile
executor_from_profile = client_from_profile


def resolve_sprkkr_profile(profile=None, *, configuration=None):
    configuration = configuration or config
    return _configured_profile(profile, configuration)
