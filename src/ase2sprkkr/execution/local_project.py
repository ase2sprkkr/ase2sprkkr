"""Reconnectable local projects use the same supervisor as SSH execution.

A project RunRecord belongs to the client. The supervisor owns only process
facts, a temporary submission ticket and a terminal outcome receipt.
"""
from __future__ import annotations

import asyncio
import os
from dataclasses import replace
from threading import Event, Lock, Thread

from .base import (
    ExecutionOutcome, ExecutionSpec, FailurePhase, OutputObservers,
    RunState, RunStatus, dispatch_cancellation,
    finish_monitor, wait_for_monitor, file_path,
)
from .file_availability import discover_local_files
from .remote_service import (
    _decode_handle, cancel_remote, follow_remote, status_remote, submit_remote,
)

_TERMINAL = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _file_path(value, directory=None):
    # Local project submission flushes file objects before publishing paths.
    return file_path(value, directory, flush=True)


def local_manifest(files, status):
    """Resolve stream availability, including partial cancelled-run output.

    The supervisor writes streams directly during execution. Cancellation and
    execution failure do not retract already written bytes.
    """
    published = (
        status.state in (RunState.FINISHED, RunState.CANCELLED) or
        (status.state == RunState.FAILED and status.failure_phase == FailurePhase.EXECUTION)
    )
    for entry in (files.stdout, files.stderr):
        if entry is None or entry.local_path is None:
            continue
        if status.failure_phase == FailurePhase.ARTIFACTS:
            state, message = "failed", status.message or "Output publication failed"
        elif published:
            exists = os.path.isfile(entry.local_path)
            state, message = ("available", None) if exists else (
                "failed", "Published output file is missing",
            )
        else:
            state, message = "skipped", None
        files = files.with_file_availability(entry.local_path, state, message)

    if status.failure_phase in (FailurePhase.LAUNCH, FailurePhase.SUBMISSION):
        return files
    return discover_local_files(files, include_streams=False)


def _local_status_from_facts(backend_id, facts, files=None):
    """Project local supervisor facts into a status without querying it again.

    Recovery already has an authoritative reply; asking the supervisor for
    status a second time risks replacing that snapshot with newer/unknown facts.
    """
    status = RunStatus(
        state=RunState(facts["state"]), backend_id=backend_id,
        return_code=facts.get("return_code"), message=facts.get("error"),
        failure_phase=facts.get("failure_phase"),
    )
    if files is not None and status.state in _TERMINAL:
        status = replace(status, files=local_manifest(files, status))
    return status


def local_status(backend_id, files=None):
    return _local_status_from_facts(
        backend_id, status_remote(backend_id), files,
    )


class _ObserverSink:
    """Adapt the supervisor's byte stream to the existing listener contract."""

    def __init__(self, observers, kind):
        self.observers = observers
        self.kind = kind

    def write(self, data):
        if data:
            self.observers.emit(data.decode("utf8", errors="replace"), self.kind)

    def flush(self):
        pass


class LocalProjectRunner:
    """Project-associated process with supervisor-based follow and reattach."""

    def __init__(
        self, executor, spec: ExecutionSpec, *, run_id=None, callback=None,
        output_listener=None, start=False, backend_id=None, state=RunState.NEW,
        failure_phase=None, failure_message=None,
    ):
        self.executor = executor
        self.spec = spec
        self.callback = callback
        self._run_id = run_id
        self._backend_id = backend_id
        self.result = ExecutionOutcome(
            output_file=_file_path(spec.output_file, spec.directory),
            input_file=_file_path(spec.input_file, spec.directory),
            error_file=_file_path(spec.error_file, spec.directory),
        )
        self._state = RunState(state)
        self._failure_phase = failure_phase
        self._failure_message = failure_message
        self._return_code = None
        self._files = None
        self._output_observers = OutputObservers(output_listener)
        self._lock = Lock()
        self._done = Event()
        self._thread = None
        self._exception = None
        self._callback_exception = None
        self._cancel_future = None
        if start:
            self.start()

    def set_manifest(self, files):
        self._files = files

    @property
    def files(self):
        return self._files

    @property
    def backend_id(self):
        return self._backend_id

    def _refresh(self):
        # A terminal result is immutable. After acknowledgment the receipt
        # may disappear while the supervisor process is still shutting down.
        if self._backend_id is None or self._state in _TERMINAL:
            return self._state
        status = local_status(self._backend_id, self._files)
        self._state = status.state
        self._failure_phase = status.failure_phase
        self._failure_message = status.message
        self._return_code = status.return_code
        self.result.return_code = self._return_code
        if status.files is not None:
            self._files = status.files
        return self._state

    @property
    def state(self):
        return self._refresh() if self._backend_id is not None else self._state

    @property
    def exception(self):
        return self._exception

    @property
    def return_code(self):
        return self._return_code

    @property
    def failure_phase(self):
        return self._failure_phase

    @property
    def failure_message(self):
        return self._failure_message

    def add_output_listener(self, callback):
        self._output_observers.add(callback)

    def add_output_error_listener(self, callback):
        self._output_observers.add_error_listener(callback)

    def _monitor(self):
        try:
            follow_remote(
                self._backend_id,
                _ObserverSink(self._output_observers, "out"),
                _ObserverSink(self._output_observers, "err"),
                stdout_path=self.result.output_file,
                stderr_path=self.result.error_file,
            )
        except (OSError, ConnectionError, RuntimeError, ValueError) as exc:
            # The follow channel is observational. A supervisor that has
            # already exited may have no socket, or no replayable spool.
            # Backend outcome is authoritative, never a guessed exit code.
            status = self._refresh()
            if status not in _TERMINAL:
                self._state = RunState.UNKNOWN
                self._failure_message = str(exc)
                return self.result
        self._refresh()
        return self.result

    def _close_owned_files(self):
        seen = set()
        for value in (self.spec.output_file, self.spec.error_file):
            if value is None or id(value) in seen:
                continue
            seen.add(id(value))
            if hasattr(value, "close"):
                value.close()

    def _thread_main(self):
        try:
            self._monitor()
        except BaseException as exc:
            self._exception = exc
            self._failure_message = str(exc)
            self._failure_phase = self._failure_phase or FailurePhase.EXECUTION
            self._state = RunState.FAILED
        finally:
            try:
                self._close_owned_files()
            except BaseException as exc:
                if self._exception is None:
                    self._exception = exc
                    self._failure_message = str(exc)
                    self._failure_phase = self._failure_phase or FailurePhase.EXECUTION
                    self._state = RunState.FAILED
            finally:
                finish_monitor(self)

    def _start_monitor(self):
        with self._lock:
            if self._thread is not None:
                return
            if self._backend_id is None:
                raise RuntimeError("Local run has no backend identity")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def _request(self):
        return {
            "executor": "local",
            "id": self._run_id,
            "detached": True,
            "directory": os.fspath(self.spec.directory),
            **self.spec.command_request(),
            "input_file": self.result.input_file,
            "stdin": _file_path(self.spec.stdin, self.spec.directory),
            "output_file": self.result.output_file,
            "error_file": self.result.error_file,
            "output_files": [
                path for value in self.spec.output_files
                if (path := _file_path(value, self.spec.directory)) is not None
            ],
            "input_files": [
                path for value in self.spec.input_files
                if (path := _file_path(value, self.spec.directory)) is not None
            ],
            "cache_files": [
                path for value in self.spec.cache_files
                if (path := _file_path(value, self.spec.directory)) is not None
            ],
        }

    def submit(self):
        if self._state != RunState.NEW:
            raise RuntimeError(f"Run is already {self._state.value}")
        if not self._run_id:
            raise RuntimeError("Project run has no run id")
        from .base import SubmissionNotStarted
        try:
            reply = submit_remote(self._request())
        except SubmissionNotStarted as exc:
            self._state = RunState.FAILED
            self._failure_phase = FailurePhase.SUBMISSION
            self._failure_message = str(exc)
            raise
        self._backend_id = reply["remote_id"]
        self._state = RunState(reply["state"])
        self._start_monitor()

    def start(self):
        if self._backend_id is None:
            return self.submit()
        return self._start_monitor()

    def wait(self):
        if self._thread is None:
            if self._state == RunState.NEW:
                raise RuntimeError("Run has not been submitted")
            raise RuntimeError("Attached run has not been started")
        return wait_for_monitor(self)

    def run(self):
        if self._thread is None:
            self.start()
        return self.wait()

    async def run_async(self):
        return await asyncio.to_thread(self.run)

    def cancel(self):
        if self._backend_id is None:
            raise RuntimeError("Run has not been submitted")
        with self._lock:
            if self._cancel_future is None:
                self._cancel_future = dispatch_cancellation(
                    lambda: cancel_remote(self._backend_id)
                )
            return self._cancel_future


def attach_local(executor, run, *, callback=None, output_listener=None, start=True):
    if run.execution.endpoint.executor != "local":
        raise ValueError("Run does not belong to the local backend")
    backend_id = run.execution.backend_id
    if _decode_handle(backend_id).get("kind") != "process":
        raise ValueError("Local project does not have a supervisor handle")
    files = run.files
    if files.stdout is None and files.stderr is None:
        raise ValueError("This run has no replayable stdout or stderr")
    spec = ExecutionSpec(
        argv=(), input_file=files.input_file,
        output_file=files.output_file,
        error_file=None if files.stderr is None else files.stderr.local_path,
        directory=run.execution.directory.local,
    )
    status = local_status(backend_id, files)
    runner = LocalProjectRunner(
        executor, spec, callback=callback, output_listener=output_listener,
        backend_id=backend_id, run_id=run.id,
        state=status.state, failure_phase=status.failure_phase,
        failure_message=status.message,
    )
    runner.set_manifest(status.files or files)
    runner._return_code = status.return_code
    runner.result.return_code = status.return_code
    if start:
        runner.start()
    return runner
