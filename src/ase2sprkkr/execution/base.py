"""Backend-independent calculation execution contracts."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional, Protocol, Sequence


class RunState(str, Enum):
    """Backend-independent state of a calculation run."""

    NEW = "new"
    QUEUED = "queued"
    RUNNING = "running"
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
class ExecutionSpec:
    """A concrete command ready to be handed to an execution backend."""

    argv: Sequence[str]
    input_file: Any
    output_file: Any
    directory: str
    stdin: Any = None
    error_file: Any = None
    input_files: Sequence[Any] = ()
    output_files: Sequence[Any] = ()


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

    def resolve_mpi(self, mpi=None):
        """Resolve an MPI setting for the machine where this executor runs.

        Local executors use the normal ASE2SPRKKR MPI autodetection.  Executors
        targeting another machine can override this method and perform the
        resolution on that target instead.
        """
        from ..configuration import mpi_runner

        return mpi_runner(mpi)

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Runner:
        raise NotImplementedError

    def attach(self, reference: RunReference, reader, callback=None) -> Runner:
        """Reconnect to a persistent run."""
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )
