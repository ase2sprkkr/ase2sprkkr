"""Backend-independent calculation execution contracts."""

from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional, Protocol, Sequence, Tuple


class MpiNotFoundWarning(UserWarning):
    """MPI autodetection was requested, but no MPI launcher was found."""


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
    UNKNOWN = "unknown"


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
    number of MPI ranks. ``mpi_warning`` controls whether failed MPI
    autodetection emits :class:`MpiNotFoundWarning`.

    The remaining fields describe remote filesystem mapping and synchronization
    for this run. They intentionally belong to the run rather than to an
    executor so executors can be shared by profiles that use the same endpoint
    but different launch/staging parameters.
    """

    mpi: Any = False
    process_count: Optional[int] = None
    mpi_warning: bool = True
    local_root: Optional[str] = None
    remote_root: Optional[str] = None
    input_sync: str = "changed"
    output_sync: str = "changed"

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


class Executor:
    """Factory/backend for calculation runners."""

    name = None

    def submission_directory(self, spec: ExecutionSpec) -> str:
        """Return this run's directory in the backend filesystem namespace."""
        return os.path.abspath(os.fspath(spec.directory))

    def resolve_mpi(self, mpi=None):
        """Convert profile MPI syntax to portable launch intent.

        No host probing is done here. The executor that actually starts the
        calculation decides how to materialize the returned launch request.
        """
        return LaunchSpec.from_mpi(mpi)

    def status(self, *references: RunReference) -> Tuple[RunStatus, ...]:
        """Return ordered status snapshots for persistent run references.

        Backends override this method directly. There is deliberately one API
        for one or many jobs so scheduler-backed executors can batch queries.
        """
        if not references:
            return ()
        raise NotImplementedError(
            f"{type(self).__name__} cannot query persistent run status"
        )

    def recover_submission(self, run_id, submission) -> Optional[RunStatus]:
        """Recover a submit whose backend reference was not persisted.

        Return a status with a reconstructed backend reference when the submit
        is found, ``None`` only when the backend can positively establish that
        it does not exist, and raise :class:`NotImplementedError` when this
        backend cannot perform such recovery. Transport/backend failures should
        propagate rather than being mistaken for absence.
        """
        raise NotImplementedError(
            f"{type(self).__name__} cannot recover an uncertain submission"
        )

    def create_persistent_run(
        self, reader, spec: ExecutionSpec, *, run_id: str, callback=None
    ) -> Runner:
        """Create a run with an ASE2SPRKKR identity known before submission.

        Backends that can expose the identity to an external scheduler or helper
        override this method. Other backends use the normal runner unchanged.
        """
        return self.create_run(reader, spec, callback=callback)

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Runner:
        raise NotImplementedError

    def attach(self, reference: RunReference, reader, callback=None) -> Runner:
        """Reconnect to a persistent run."""
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )
