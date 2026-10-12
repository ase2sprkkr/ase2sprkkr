"""Backend-independent calculation execution contracts."""

from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional, Protocol, Sequence, Tuple


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
    cache_files: Sequence[Any] = ()


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
        """Resolve an MPI setting for the machine where this executor runs.

        ``None`` means that the selected execution profile did not request MPI;
        it must not fall back to another global configuration authority.
        Executors targeting another machine can override this method and resolve
        non-trivial MPI settings on that target instead.
        """
        if mpi is None:
            return None

        from ..configuration import mpi_runner

        return mpi_runner(mpi)

    def create_run(self, reader, spec: ExecutionSpec, callback=None) -> Runner:
        raise NotImplementedError

    def attach(self, reference: RunReference, reader, callback=None) -> Runner:
        """Reconnect to a persistent run."""
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )
