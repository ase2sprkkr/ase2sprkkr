"""Backend-independent calculation execution contracts."""

from __future__ import annotations

import os
from dataclasses import dataclass, replace
from enum import Enum
from typing import Any, Optional, Protocol, Sequence, Tuple


class MPIRunnerNotFoundWarning(UserWarning):
    """Warning emitted when an explicitly requested MPI launcher is unavailable."""


class RunState(str, Enum):
    """Backend-independent state of a calculation run."""

    NEW = "new"
    SUBMITTING = "submitting"
    QUEUED = "queued"
    RUNNING = "running"
    FINISHING = "finishing"
    FINISHED = "finished"
    FAILED = "failed"
    CANCELLED = "cancelled"
    LOST = "lost"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class RunTarget:
    """Filesystem resources owned by a run on one execution target.

    ``files`` contains paths that the run may modify. Two targets conflict only
    when they refer to the same filesystem namespace and own at least one of
    the same files.
    """

    location: str
    files: Tuple[str, ...] = ()

    def __post_init__(self):
        object.__setattr__(self, "files", tuple(sorted(set(self.files))))

    def conflicts_with(self, other: "RunTarget") -> bool:
        return self.location == other.location and bool(
            set(self.files) & set(other.files)
        )


@dataclass(frozen=True)
class RunReference:
    """Serializable information needed to reconnect to a persistent run."""

    executor: str
    native_id: str
    computer: Optional[str] = None
    directory: Optional[str] = None
    input_file: Optional[str] = None
    output_file: Optional[str] = None
    error_file: Optional[str] = None
    target: Optional[RunTarget] = None

    def __post_init__(self):
        if isinstance(self.target, dict):
            object.__setattr__(self, "target", RunTarget(**self.target))


@dataclass(frozen=True)
class RunStatus:
    """One backend status snapshot for a persistent run reference."""

    reference: RunReference
    state: RunState
    return_code: Optional[int] = None
    message: Optional[str] = None


@dataclass(frozen=True)
class LaunchSpec:
    """Portable launch and staging parameters for one concrete run.

    ``mpi`` is ``False`` for serial execution, ``True`` when MPI is required,
    ``"auto"`` when a local executor may fall back to serial execution, or an
    explicit launcher argv tuple. ``process_count`` is an optional requested
    number of MPI ranks. ``mpi_warning`` controls whether failure to detect an
    explicitly requested default MPI launcher emits a warning.

    The remaining fields describe remote filesystem mapping and synchronization
    for this run. They intentionally belong to the run rather than to an
    executor so executors can be shared by profiles that use the same endpoint
    but different launch/staging parameters.
    """

    mpi: Any = False
    process_count: Optional[int] = None
    local_root: Optional[str] = None
    remote_root: Optional[str] = None
    input_sync: str = "changed"
    output_sync: str = "changed"
    mpi_warning: bool = True

    def __post_init__(self):
        if isinstance(self.mpi, list):
            object.__setattr__(self, "mpi", tuple(self.mpi))
        if self.process_count is not None and self.process_count < 1:
            raise ValueError("process_count must be positive")
        for name in ("local_root", "remote_root"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, os.fspath(value))
        for name in ("input_sync", "output_sync"):
            value = getattr(self, name)
            object.__setattr__(self, name, getattr(value, "value", value))

    def __bool__(self):
        """Whether this run requests MPI."""
        return self.mpi is not False and self.mpi is not None

    @classmethod
    def from_mpi(cls, mpi):
        """Translate the existing ``mpi`` profile syntax without probing a host."""
        if isinstance(mpi, cls):
            return mpi
        if mpi is None or mpi is False:
            return cls()
        if isinstance(mpi, bool):
            return cls(mpi=mpi)
        if isinstance(mpi, int):
            if mpi < 1:
                raise ValueError("MPI process count must be positive")
            return cls(mpi=True, process_count=mpi)
        if isinstance(mpi, str):
            return cls(mpi="auto" if mpi == "auto" else (mpi,))
        return cls(mpi=tuple(mpi))


def _pathlike(value):
    """Normalize path-like values without touching open file objects."""
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


@dataclass(frozen=True)
class ExecutionSpec:
    """A command and launch intent ready to be handed to an execution backend."""

    argv: Sequence[str]
    input_file: Any
    output_file: Any
    directory: str
    stdin: Any = None
    error_file: Any = None
    input_files: Sequence[Any] = ()
    output_files: Sequence[Any] = ()
    cache_files: Sequence[Any] = ()
    mpi_argv: Optional[Sequence[str]] = None
    launch: LaunchSpec = LaunchSpec()

    def __post_init__(self):
        object.__setattr__(self, "directory", os.fspath(self.directory))
        for name in ("input_file", "output_file", "stdin", "error_file"):
            object.__setattr__(self, name, _pathlike(getattr(self, name)))
        for name in ("input_files", "output_files", "cache_files"):
            object.__setattr__(
                self,
                name,
                tuple(_pathlike(value) for value in getattr(self, name)),
            )
        object.__setattr__(self, "argv", tuple(os.fspath(value) for value in self.argv))
        if self.mpi_argv is not None:
            object.__setattr__(
                self,
                "mpi_argv",
                tuple(os.fspath(value) for value in self.mpi_argv),
            )


class Runner(Protocol):
    """One concrete calculation execution."""

    result: Any

    @property
    def state(self) -> RunState: ...

    @property
    def reference(self) -> Optional[RunReference]: ...

    def submit(self) -> None: ...

    def wait(self) -> Any: ...

    def run(self) -> Any: ...

    def cancel(self) -> None: ...

    def add_output_listener(self, callback) -> None: ...


def _resource_path(value, directory):
    if value is None:
        return None
    value = getattr(value, "name", value)
    if isinstance(value, int):
        # Anonymous temporary files and subprocess special values do not have a
        # filesystem name that can collide with another calculation.
        return None
    try:
        path = os.fspath(value)
    except TypeError:
        return None
    if not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


class Executor:
    """Factory/backend for calculation runners."""

    name = None

    def execution_target(self, spec: ExecutionSpec) -> RunTarget:
        """Return the concrete resources that can collide with another run."""
        directory = os.path.abspath(spec.directory)
        files = set()
        for value in (
            spec.output_file,
            spec.error_file,
            *spec.output_files,
            *spec.cache_files,
        ):
            path = _resource_path(value, directory)
            if path is not None:
                files.add(path)
        return RunTarget("local", tuple(files))

    def resolve_mpi(self, mpi=None):
        """Convert profile MPI syntax to portable launch intent.

        No host probing is done here. The executor that actually starts the
        calculation decides how to materialize the returned launch request.
        """
        return LaunchSpec.from_mpi(mpi)

    def status(self, *references: RunReference) -> Tuple[RunStatus, ...]:
        """Return status snapshots for all supplied persistent references.

        There is deliberately one API for one or many jobs. Backends may batch
        the underlying query, for example into one scheduler status command.
        The returned tuple always follows the input order.
        """
        from .status import status

        return status(self, *references)

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Runner:
        raise NotImplementedError

    def attach(self, reference: RunReference, reader, callback=None) -> Runner:
        """Reconnect to a persistent run."""
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )
