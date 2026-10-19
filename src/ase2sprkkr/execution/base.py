"""Backend-independent calculation execution contracts."""

from __future__ import annotations

import asyncio
import os
import warnings
from dataclasses import dataclass
from enum import Enum
from threading import Event, Lock
from typing import Any, Optional, Protocol, Sequence, Tuple

from ..common.decorators import class_property


class MpiNotFoundWarning(UserWarning):
    """MPI autodetection was requested, but no MPI launcher was found."""


class OutputObserverWarning(RuntimeWarning):
    """An output observer failed while the calculation itself keeps running."""


class SubmissionNotStarted(RuntimeError):
    """Submission failed before the backend could possibly accept the run."""


class OutputObservers:
    """Deliver output without allowing observers to control execution outcome."""

    def __init__(self, listener=None):
        self._listeners = [] if listener is None else [listener]
        self._error_listeners = []
        self._lock = Lock()

    def add(self, callback):
        with self._lock:
            self._listeners.append(callback)

    def add_error_listener(self, callback):
        with self._lock:
            self._error_listeners.append(callback)

    def emit(self, text, kind):
        with self._lock:
            listeners = tuple(self._listeners)
        for callback in listeners:
            try:
                callback(text, kind)
            except Exception as exc:
                with self._lock:
                    try:
                        self._listeners.remove(callback)
                    except ValueError:
                        pass
                    handlers = tuple(self._error_listeners)
                if not handlers:
                    warnings.warn(
                        f"Output observer {callback!r} failed: {exc}",
                        OutputObserverWarning,
                        stacklevel=4,
                    )
                    continue
                for handler in handlers:
                    try:
                        handler(callback, exc)
                    except Exception as handler_exc:
                        warnings.warn(
                            f"Output observer error handler {handler!r} failed: "
                            f"{handler_exc}",
                            OutputObserverWarning,
                            stacklevel=4,
                        )


class RunState(str, Enum):
    """Backend-independent state of a calculation run."""

    NEW = "new"
    SUBMITTING = "submitting"
    QUEUED = "queued"
    RUNNING = "running"
    FINALIZING = "finalizing"
    FINISHING = "finalizing"
    FINISHED = "finished"
    EXECUTION_FAILED = "execution_failed"
    FAILED = "execution_failed"
    FINALIZATION_FAILED = "finalization_failed"
    CANCELLED = "cancelled"
    UNKNOWN = "unknown"


class FailurePhase(str, Enum):
    """Stage in which a terminal failure occurred.

    This is diagnostic lifecycle detail, not a second state machine. In
    particular, ``LAUNCH`` means that no calculation process ever reached the
    running state, so no output parsing/finalization belongs to that run.
    """

    SUBMISSION = "submission"
    LAUNCH = "launch"
    EXECUTION = "execution"
    FINALIZATION = "finalization"


class ExecutionResources:
    """Closable resources whose ownership is transferred to a calculation run.

    Closing is idempotent and thread-safe. Duplicate objects are closed only
    once, which matters when two execution roles reference the same file handle.
    """

    def __init__(self, *resources):
        self._resources = tuple(resources)
        self._closed = False
        self._lock = Lock()

    @property
    def closed(self):
        with self._lock:
            return self._closed

    def close(self, *, suppress_errors=False):
        with self._lock:
            if self._closed:
                return
            self._closed = True
            resources = self._resources
            self._resources = ()

        first_error = None
        seen = set()
        for resource in resources:
            if resource is None or id(resource) in seen:
                continue
            seen.add(id(resource))
            close = getattr(resource, "close", None)
            if close is None or getattr(resource, "closed", False):
                continue
            try:
                close()
            except BaseException as exc:
                if first_error is None:
                    first_error = exc

        if first_error is not None and not suppress_errors:
            raise first_error


class ResourceOwnedRunner:
    """Runner proxy owning resources and the public submit/cancel boundary.

    Resources remain owned by the caller until this proxy is constructed. From
    that point they are released on completion, on a synchronous submission
    failure, or explicitly through :meth:`close` while the run is still new.

    Cancellation is latched and thread-safe even before submission. If
    :meth:`cancel` wins the race with :meth:`submit`, the backend runner is not
    submitted and ``submit()`` raises :class:`asyncio.CancelledError`. If the
    request arrives while backend submission is already in progress, cancellation
    is dispatched exactly once after that submission boundary completes.
    """

    def __init__(self, runner, resources: ExecutionResources):
        self._runner = runner
        self._resources = resources
        self._submission_lock = Lock()
        self._cancel_requested = Event()
        # A proxy may also be constructed around an already-started runner.
        # Only NEW unambiguously means the backend submission boundary has not
        # been crossed; for every other state cancellation must be delegated.
        self._submitted = runner.state != RunState.NEW
        self._cancel_dispatched = False

    def __getattr__(self, name):
        return getattr(self._runner, name)

    @property
    def result(self):
        return self._runner.result

    @property
    def state(self):
        if self._cancel_requested.is_set():
            with self._submission_lock:
                if not self._submitted:
                    return RunState.CANCELLED
        return self._runner.state

    @property
    def reference(self):
        return self._runner.reference

    @property
    def exception(self):
        return self._runner.exception

    def _release_after_failure(self):
        self._resources.close(suppress_errors=True)

    def _raise_if_cancelled_before_submit(self):
        if not self._cancel_requested.is_set():
            return
        with self._submission_lock:
            if not self._submitted:
                raise asyncio.CancelledError("Run was cancelled before submission")

    def submit(self):
        with self._submission_lock:
            if self._cancel_requested.is_set():
                self._resources.close()
                raise asyncio.CancelledError("Run was cancelled before submission")
            try:
                result = self._runner.submit()
            except BaseException:
                self._release_after_failure()
                raise
            self._submitted = True
            cancel_after_submit = self._cancel_requested.is_set()

        if cancel_after_submit:
            self.cancel()
        return result

    def wait(self):
        self._raise_if_cancelled_before_submit()
        try:
            result = self._runner.wait()
        except BaseException:
            self._release_after_failure()
            raise
        self._resources.close()
        return result

    def run(self):
        self._raise_if_cancelled_before_submit()
        if self.state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        self._raise_if_cancelled_before_submit()
        if self.state == RunState.NEW:
            self.submit()
        try:
            run_async = getattr(self._runner, "run_async", None)
            if run_async is not None:
                result = await run_async()
            else:
                result = await asyncio.to_thread(self._runner.wait)
        except BaseException:
            self._release_after_failure()
            raise
        self._resources.close()
        return result

    def cancel(self):
        self._cancel_requested.set()
        with self._submission_lock:
            if not self._submitted:
                self._resources.close()
                return None
            if self._cancel_dispatched:
                return None
            result = self._runner.cancel()
            self._cancel_dispatched = True
            return result

    def close(self):
        """Release resources belonging to a run that has not been started."""
        if self._resources.closed:
            return
        state = self.state
        if state != RunState.NEW:
            raise RuntimeError(
                f"Cannot close execution resources while run is {state.value}"
            )
        self._resources.close()

    def add_output_listener(self, callback):
        return self._runner.add_output_listener(callback)

    def add_output_error_listener(self, callback):
        return self._runner.add_output_error_listener(callback)


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
    failure_phase: Optional[FailurePhase] = None

    def __post_init__(self):
        if not isinstance(self.state, RunState):
            object.__setattr__(self, "state", RunState(self.state))
        if self.failure_phase is not None and not isinstance(
            self.failure_phase, FailurePhase
        ):
            object.__setattr__(
                self, "failure_phase", FailurePhase(self.failure_phase)
            )


@dataclass(frozen=True)
class MpiSpec:
    """Portable MPI launch semantics independent of the execution host.

    ``autodetect`` requests host-specific launcher discovery. ``required``
    distinguishes mandatory MPI from ``auto`` mode, which may fall back to
    serial execution. An explicit ``launcher`` disables autodetection.
    ``process_count`` is kept portable until the launcher is materialized.
    """

    required: bool = False
    autodetect: bool = False
    launcher: Tuple[str, ...] = ()
    process_count: Optional[int] = None
    warning: bool = True

    def __post_init__(self):
        object.__setattr__(
            self, "launcher", tuple(os.fspath(part) for part in self.launcher)
        )
        if self.process_count is not None and self.process_count < 1:
            raise ValueError("process_count must be positive")
        if self.launcher and self.autodetect:
            raise ValueError("explicit MPI launcher cannot be combined with autodetection")
        if self.launcher and not self.required:
            raise ValueError("explicit MPI launcher must be required")
        if self.required and not self.autodetect and not self.launcher:
            raise ValueError("required MPI needs autodetection or an explicit launcher")
        if not self.requested and self.process_count is not None:
            raise ValueError("process_count requires MPI")

    @classmethod
    def from_value(cls, mpi, *, process_count=None, warning=True):
        """Translate public ``mpi`` syntax into normalized portable semantics."""
        if isinstance(mpi, cls):
            if process_count is not None or warning is not True:
                raise ValueError("MPI options cannot override an existing MpiSpec")
            return mpi
        if mpi is None or mpi is False:
            if process_count is not None:
                raise ValueError("process_count requires MPI")
            return cls(warning=warning)
        if isinstance(mpi, bool):
            return cls(
                required=True,
                autodetect=True,
                process_count=process_count,
                warning=warning,
            )
        if isinstance(mpi, int):
            if process_count is not None:
                raise ValueError("MPI process count was specified twice")
            return cls(
                required=True,
                autodetect=True,
                process_count=mpi,
                warning=warning,
            )
        if mpi == "auto":
            return cls(
                autodetect=True,
                process_count=process_count,
                warning=warning,
            )
        if isinstance(mpi, str):
            launcher = (mpi,)
        else:
            launcher = tuple(os.fspath(part) for part in mpi)
        return cls(
            required=True,
            launcher=launcher,
            process_count=process_count,
            warning=warning,
        )

    @property
    def requested(self) -> bool:
        """Whether the calculation requests an MPI launch."""
        return self.required or self.autodetect or bool(self.launcher)

    @property
    def allows_serial_fallback(self) -> bool:
        """Whether failure to find MPI may fall back to serial execution."""
        return self.autodetect and not self.required


@dataclass(frozen=True)
class LaunchSpec:
    """Portable launch and staging parameters for one concrete run.

    MPI semantics are represented by :class:`MpiSpec`; host-independent MPI
    intent therefore stays grouped instead of being spread over parallel launch
    fields. The remaining fields describe remote filesystem mapping and
    synchronization for this run.
    """

    mpi: MpiSpec = MpiSpec()
    local_root: Optional[str] = None
    remote_root: Optional[str] = None
    input_sync: str = "changed"
    output_sync: str = "changed"

    def __post_init__(self):
        mpi = self.mpi
        if isinstance(mpi, dict):
            mpi = MpiSpec(**mpi)
        elif not isinstance(mpi, MpiSpec):
            raise TypeError("LaunchSpec.mpi must be an MpiSpec")
        object.__setattr__(self, "mpi", mpi)

        for name in ("local_root", "remote_root"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, os.fspath(value))
        for name in ("input_sync", "output_sync"):
            value = getattr(self, name)
            object.__setattr__(self, name, getattr(value, "value", value))

    def __bool__(self):
        """Whether this run requests MPI."""
        return self.mpi.requested

    @classmethod
    def from_mpi(cls, mpi):
        """Translate public ``mpi`` syntax without probing an execution host."""
        if isinstance(mpi, cls):
            return mpi
        return cls(mpi=MpiSpec.from_value(mpi))


def _pathlike(value):
    """Normalize path-like values without touching open file objects."""
    from ..input_parameters.task_definition import TaskFile
    if isinstance(value, TaskFile):
        return value
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

    @property
    def exception(self) -> Optional[BaseException]: ...

    def submit(self) -> None: ...

    def wait(self) -> Any: ...

    def run(self) -> Any: ...

    def cancel(self) -> None: ...

    def add_output_listener(self, callback) -> None: ...

    def add_output_error_listener(self, callback) -> None: ...


class Executor:
    """Factory/backend for calculation runners."""

    name = None

    @class_property
    def supports_result_file_retrieval(cls):
        """Whether this executor can restore files belonging to a stored run."""
        return cls.fetch_result_files is not Executor.fetch_result_files

    def fetch_result_files(self, reference: RunReference, files):
        """Restore requested local files belonging to a persistent run.

        ``files`` contains the concrete historical local paths, or pre-capture
        patterns, selected by the persistence layer. Backends only map and
        retrieve those requested files; they must not rediscover a historical
        result by re-expanding their own current output manifests.
        """
        raise NotImplementedError(
            f"{type(self).__name__} cannot retrieve persistent result files"
        )

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

    def recover_submission(
        self, run_id, submission, *, submitted_at=None
    ) -> Optional[RunStatus]:
        """Recover a submit whose backend reference was not persisted.

        ``submitted_at`` is the timestamp at which ASE2SPRKKR reserved the
        submission. Scheduler backends may use it to bound accounting-history
        queries after a job has left the live queue.

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

    def attach(
        self,
        reference: RunReference,
        reader,
        callback=None,
        output_listener=None,
        start=True,
    ) -> Runner:
        """Reconnect to a persistent run.

        ``output_listener`` is part of runner construction, while ``start``
        controls the lifecycle side effect. Backends must install the listener
        before honoring ``start=True`` so immediately replayed output is never
        lost.
        """
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )