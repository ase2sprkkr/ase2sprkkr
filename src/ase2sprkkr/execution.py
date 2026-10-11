"""Backend-independent execution primitives.

This module deliberately knows nothing about SPR-KKR executable naming, MPI
launchers or input syntax.  Those details belong to :mod:`sprkkr.profile`.
"""

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
    directory: Optional[str] = None
    input_file: Optional[str] = None
    output_file: Optional[str] = None


@dataclass(frozen=True)
class ExecutionSpec:
    """A concrete command ready to be handed to an execution backend."""

    argv: Sequence[str]
    input_file: Any
    output_file: Any
    directory: str
    stdin: Any = None


class Run(Protocol):
    """A runnable or running calculation."""

    result: Any

    @property
    def state(self) -> RunState: ...

    @property
    def reference(self) -> Optional[RunReference]: ...

    def run(self) -> Any: ...

    def cancel(self) -> None: ...


class Executor:
    """Base class for calculation execution backends."""

    name = None

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Run:
        raise NotImplementedError

    def attach(self, reference: RunReference, **kwargs) -> Run:
        """Reconnect to a persistent run.

        The local executor cannot do this.  Scheduler-backed executors can
        implement it using the native scheduler id stored in ``reference``.
        """
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )


class LocalExecutor(Executor):
    """Execute an already prepared command on the current machine."""

    name = "local"

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Run:
        return reader.create_process(
            list(spec.argv),
            spec.output_file,
            stdin=spec.stdin,
            input_file=spec.input_file.name,
            callback=callback,
        )
