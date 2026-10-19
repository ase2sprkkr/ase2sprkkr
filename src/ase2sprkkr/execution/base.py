"""Backend-independent calculation execution contracts."""

from __future__ import annotations

import asyncio
import os
import warnings
from concurrent.futures import Future
from dataclasses import dataclass
from enum import Enum
from threading import Event, Lock, Thread
from typing import Any, Optional, Protocol, Sequence, Tuple

from ..common.decorators import class_property
from .run_model import Execution, FileManifest


def file_path(value, directory=None, *, flush=False, required=False):
    """Resolve a named file or path, returning None for unnamed descriptors."""
    if value is None:
        return None
    if flush and hasattr(value, "flush"):
        value.flush()
    value = getattr(value, "name", value)
    if isinstance(value, int):
        if required:
            raise ValueError("A named file or path is required")
        return None
    try:
        path = os.fspath(value)
    except TypeError as exc:
        if required:
            raise ValueError("A named file or path is required") from exc
        return None
    if directory is not None and not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


def dispatch_cancellation(action):
    """Run a blocking backend cancellation off the calling/GUI thread.

    The returned Future exposes errors, while the normal monitor remains
    authoritative for the terminal run state.
    """
    future = Future()

    def worker():
        if future.set_running_or_notify_cancel():
            try:
                future.set_result(action())
            except BaseException as exc:
                future.set_exception(exc)

    Thread(target=worker, daemon=True).start()
    return future


def start_monitor(runner, *, missing_backend=None):
    """Start a runner's monitor thread only once, holding its lifecycle lock."""
    with runner._lock:
        if runner._thread is not None:
            return
        if missing_backend is not None and runner._backend_id is None:
            raise RuntimeError(missing_backend)
        runner._thread = Thread(target=runner._thread_main, daemon=True)
        runner._thread.start()


def request_cancellation(runner, action):
    """Dispatch one cancellation per runner and return the same Future on retries."""
    with runner._lock:
        if runner._cancel_future is None:
            runner._cancel_future = dispatch_cancellation(action)
        return runner._cancel_future


def finish_monitor(runner):
    """Deliver the callback and unblock waiters after backend-specific cleanup."""
    try:
        if runner.callback is not None:
            runner.callback(runner.result)
    except BaseException as exc:
        runner._callback_exception = exc
    finally:
        runner._done.set()


def wait_for_monitor(runner):
    """Propagate execution errors before callback errors after completion."""
    runner._done.wait()
    if runner._exception is not None:
        raise runner._exception
    if runner._callback_exception is not None:
        raise runner._callback_exception
    return runner.result


class MpiNotFoundWarning(UserWarning):
    """MPI autodetection was requested, but no MPI launcher was found."""


class OutputObserverWarning(RuntimeWarning):
    """An output observer failed while the calculation itself keeps running."""


class SubmissionNotStarted(RuntimeError):
    """Submission failed before the backend could possibly accept the run."""


class OutputObservers:
    """Deliver raw output without allowing observers to affect execution."""

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
    """Single authoritative lifecycle state of one calculation run."""

    NEW = "new"
    SUBMITTING = "submitting"
    QUEUED = "queued"
    RUNNING = "running"
    DOWNLOADING = "downloading"
    POSTPROCESSING = "postprocessing"
    FINISHED = "finished"
    FAILED = "failed"
    CANCELLED = "cancelled"
    UNKNOWN = "unknown"


class FailurePhase(str, Enum):
    """Phase in which a run entered :attr:`RunState.FAILED`."""

    SUBMISSION = "submission"
    LAUNCH = "launch"
    EXECUTION = "execution"
    ARTIFACTS = "artifacts"
    DOWNLOAD = "download"
    POSTPROCESSING = "postprocessing"


@dataclass
class ExecutionOutcome:
    """Raw execution facts available before task-specific interpretation."""

    output_file: Any = None
    input_file: Any = None
    error_file: Any = None
    return_code: Optional[int] = None


class ExecutionResources:
    """Closable resources whose ownership is transferred to a calculation run."""

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
    """Runner proxy owning resources and the public submit/cancel boundary."""

    def __init__(self, runner, resources: ExecutionResources):
        self._runner = runner
        self._resources = resources
        self._submission_lock = Lock()
        self._cancel_requested = Event()
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
    def backend_id(self):
        return self._runner.backend_id

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
class RunStatus:
    """Backend status facts; the client owns durable RunRecord state."""

    state: RunState
    backend_id: Optional[str] = None
    execution: Optional[Execution] = None
    files: Optional[FileManifest] = None
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
    """Portable MPI launch semantics independent of the execution host."""

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
        launcher = (mpi,) if isinstance(mpi, str) else tuple(os.fspath(part) for part in mpi)
        return cls(
            required=True,
            launcher=launcher,
            process_count=process_count,
            warning=warning,
        )

    @property
    def requested(self) -> bool:
        return self.required or self.autodetect or bool(self.launcher)

    @property
    def allows_serial_fallback(self) -> bool:
        return self.autodetect and not self.required


@dataclass(frozen=True)
class LaunchSpec:
    """Portable launch and staging parameters for one concrete run."""

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
        return self.mpi.requested

    @classmethod
    def from_mpi(cls, mpi):
        if isinstance(mpi, cls):
            return mpi
        return cls(mpi=MpiSpec.from_value(mpi))


def _pathlike(value):
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


@dataclass(frozen=True)
class ExecutionSpec:
    """Command, files and launch intent handed to an execution backend."""

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
    # The MPI command may have a different stdin contract from serial argv.
    # Generic executors preserve stdin unless the caller explicitly opts out.
    mpi_uses_stdin: bool = True

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
    """One concrete calculation run."""

    result: Any

    @property
    def state(self) -> RunState: ...

    @property
    def backend_id(self) -> Optional[str]: ...

    @property
    def exception(self) -> Optional[BaseException]: ...

    def submit(self) -> None: ...
    def wait(self) -> Any: ...
    def run(self) -> Any: ...
    def cancel(self) -> None: ...
    def add_output_listener(self, callback) -> None: ...
    def add_output_error_listener(self, callback) -> None: ...


class Executor:
    """Execution backend. It never knows task readers or TaskResult types."""

    name = None

    @class_property
    def supports_result_file_retrieval(cls):
        return cls.fetch_result_files is not Executor.fetch_result_files

    def fetch_result_files(self, run, files):
        raise NotImplementedError(
            f"{type(self).__name__} cannot retrieve persistent result files"
        )

    def submission_directory(self, spec: ExecutionSpec) -> str:
        return os.path.abspath(os.fspath(spec.directory))

    def resolve_mpi(self, mpi=None):
        return LaunchSpec.from_mpi(mpi)

    def status(self, *runs) -> Tuple[RunStatus, ...]:
        if not runs:
            return ()
        raise NotImplementedError(
            f"{type(self).__name__} cannot query persistent run status"
        )

    def recover_submission(
        self, run_id, submission, *, submitted_at=None
    ) -> Optional[RunStatus]:
        raise NotImplementedError(
            f"{type(self).__name__} cannot recover an uncertain submission"
        )

    def create_project_run(
        self, spec: ExecutionSpec, *, run_id: str, callback=None
    ) -> Runner:
        return self.create_run(spec, callback=callback)

    def create_run(self, spec: ExecutionSpec, callback=None) -> Runner:
        raise NotImplementedError

    def attach(
        self,
        run,
        callback=None,
        output_listener=None,
        start=True,
    ) -> Runner:
        raise NotImplementedError(
            f"{type(self).__name__} cannot attach to an existing run"
        )
