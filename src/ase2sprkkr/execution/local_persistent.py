"""Reconnectable local execution using the detached execution helper."""

from __future__ import annotations

import asyncio
import json
import os
from dataclasses import asdict
from threading import Event, Lock, Thread

from .base import (
    ExecutionOutcome,
    ExecutionSpec,
    FailurePhase,
    OutputObservers,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .detached import (
    _request_path,
    cancel_detached,
    detached_failure_phase,
    detached_launch_failed,
    detached_run_state,
    detached_stream_file,
    status_detached,
    submit_detached,
)

_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _file_path(value, directory=None):
    if value is None:
        return None
    if hasattr(value, "flush"):
        value.flush()
    value = getattr(value, "name", value)
    if isinstance(value, int):
        return None
    try:
        path = os.fspath(value)
    except TypeError:
        return None
    if directory is not None and not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


def detached_local_status(run):
    execution = run.execution
    try:
        status = status_detached(execution.directory.local, execution.backend_id)
    except (FileNotFoundError, json.JSONDecodeError):
        return None
    return RunStatus(
        state=detached_run_state(status.get("state", "unknown")),
        backend_id=execution.backend_id,
        return_code=status.get("return_code"),
        message=status.get("error"),
        failure_phase=detached_failure_phase(status),
    )


def recover_local_submission(run_id, run):
    directory = run.execution.directory.local
    try:
        state = status_detached(directory, run_id)
    except FileNotFoundError:
        return None
    return RunStatus(
        state=detached_run_state(state.get("state", "unknown")),
        backend_id=run_id,
        return_code=state.get("return_code"),
        message=state.get("error"),
        failure_phase=detached_failure_phase(state),
    )


class LocalPersistentRunner:
    """Detached local execution with durable helper-owned artifact publication."""

    def __init__(
        self,
        executor,
        spec: ExecutionSpec,
        *,
        run_id=None,
        callback=None,
        output_listener=None,
        start=False,
        backend_id=None,
        state=RunState.NEW,
        failure_phase=None,
        failure_message=None,
    ):
        self.executor = executor
        self.spec = spec
        self.callback = callback
        self._run_id = run_id or backend_id
        self._backend_id = backend_id
        self.result = ExecutionOutcome(
            output_file=_file_path(spec.output_file, spec.directory),
            input_file=_file_path(spec.input_file, spec.directory),
            error_file=_file_path(spec.error_file, spec.directory),
        )

        self._state = RunState(state)
        self._return_code = None
        self._failure_phase = (
            None if failure_phase is None else FailurePhase(failure_phase)
        )
        self._failure_message = failure_message
        self._output_observers = OutputObservers(output_listener)
        self._lock = Lock()
        self._done = Event()
        self._thread = None
        self._exception = None
        self._callback_exception = None

        if start:
            self.start()

    @property
    def backend_id(self):
        return self._backend_id

    @property
    def state(self):
        if self._backend_id is None or self._state in _FINAL_STATES:
            return self._state
        from .run_model import Execution, WorkingDirectory
        from .endpoint import ExecutorEndpoint
        from types import SimpleNamespace
        run = SimpleNamespace(execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(os.path.abspath(self.spec.directory)),
            self._backend_id,
        ))
        status = detached_local_status(run)
        self._state = RunState.UNKNOWN if status is None else status.state
        if status is not None:
            self._failure_phase = status.failure_phase
            self._failure_message = status.message
            if status.return_code is not None:
                self._return_code = status.return_code
                self.result.return_code = status.return_code
        return self._state

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

    def _emit_output(self, data, kind):
        if data:
            self._output_observers.emit(data.decode("utf8"), kind)

    @staticmethod
    def _read_new(path, offset):
        if not path:
            return b""
        try:
            with open(path, "rb") as source:
                source.seek(offset)
                return source.read()
        except FileNotFoundError:
            return b""

    async def _monitor(self):
        offsets = {"out": 0, "err": 0}
        while True:
            status = status_detached(
                self.spec.directory, self._backend_id
            )
            self._failure_phase = detached_failure_phase(status)
            self._failure_message = status.get("error")
            if detached_launch_failed(status):
                self._state = RunState.FAILED
                self._return_code = None
                self._failure_phase = FailurePhase.LAUNCH
                return self.result

            state = detached_run_state(status.get("state", "unknown"))
            started = status.get("started_at") is not None
            had_data = False
            if started:
                seen = set()
                for kind, key in (("out", "output_file"), ("err", "error_file")):
                    path = detached_stream_file(status, key)
                    if not path or path in seen:
                        continue
                    seen.add(path)
                    data = self._read_new(path, offsets[kind])
                    if not data:
                        continue
                    offsets[kind] += len(data)
                    self._emit_output(data, kind)
                    had_data = True

            if state in _FINAL_STATES and not had_data:
                self._state = state
                self._return_code = status.get("return_code")
                if self._return_code is None and started:
                    self._return_code = 0 if state == RunState.FINISHED else 1
                self.result.return_code = self._return_code
                return self.result

            if not had_data:
                await asyncio.sleep(0.05)

    def _close_owned_files(self):
        closed = set()
        for file in (self.spec.output_file, self.spec.error_file):
            if file is None or id(file) in closed:
                continue
            closed.add(id(file))
            if hasattr(file, "close"):
                file.close()

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
            if self._failure_message is None:
                self._failure_message = str(exc)
            if self._failure_phase is None:
                self._failure_phase = FailurePhase.EXECUTION
            if self._state != RunState.CANCELLED:
                self._state = RunState.FAILED
        finally:
            self._close_owned_files()
            try:
                if self.callback:
                    self.callback(self.result)
            except BaseException as exc:
                self._callback_exception = exc
            self._done.set()

    def _start_monitor(self):
        with self._lock:
            if self._thread is not None:
                return
            if self._backend_id is None:
                raise RuntimeError("Run has no backend reference")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def start(self):
        if self._backend_id is None:
            return self.submit()
        return self._start_monitor()

    def _request(self):
        directory = os.path.abspath(os.fspath(self.spec.directory))
        return {
            "id": self._run_id,
            "executor": "local",
            "directory": directory,
            "argv": list(self.spec.argv),
            "mpi_argv": None
            if self.spec.mpi_argv is None
            else list(self.spec.mpi_argv),
            "launch": asdict(self.spec.launch),
            "input_file": self.result.input_file,
            "stdin": _file_path(self.spec.stdin, directory),
            # The detached helper is the single owner of durable publication.
            # It writes run-specific launch streams and atomically replaces these
            # canonical paths before publishing FINISHED.
            "output_file": self.result.output_file,
            "error_file": self.result.error_file,
        }

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")
            if not self._run_id:
                raise RuntimeError("Persistent local run has no run id")
            try:
                request = self._request()
                status = submit_detached(request)
            except OSError as exc:
                self._failure_phase = FailurePhase.SUBMISSION
                self._failure_message = str(exc)
                self._state = RunState.FAILED
                raise SubmissionNotStarted(str(exc)) from exc
            if detached_launch_failed(status):
                self._state = RunState.FAILED
                self._failure_phase = FailurePhase.LAUNCH
                self._failure_message = status.get("error")
                detail = self._failure_message or "unknown error"
                raise SubmissionNotStarted(
                    f"Detached worker launch failed: {detail}"
                )
            self._state = detached_run_state(status.get("state", "submitting"))
            self._return_code = status.get("return_code")
            self.result.return_code = self._return_code
            self._failure_phase = detached_failure_phase(status)
            self._failure_message = status.get("error")
            self._backend_id = self._run_id
        self._start_monitor()

    def wait(self):
        if self._thread is None:
            if self._state == RunState.NEW:
                raise RuntimeError("Run has not been submitted")
            raise RuntimeError("Attached run has not been started")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        if self._callback_exception is not None:
            raise self._callback_exception
        return self.result

    def run(self):
        if self._thread is None:
            self.start()
        return self.wait()

    async def run_async(self):
        return await asyncio.to_thread(self.run)

    def cancel(self):
        if self._backend_id is None:
            raise RuntimeError("Run has not been submitted")
        state = cancel_detached(
            self.spec.directory, self._backend_id
        )
        self._state = detached_run_state(state.get("state", "cancelled"))
        self._failure_phase = detached_failure_phase(state)
        self._failure_message = state.get("error")
        if state.get("return_code") is not None:
            self._return_code = state["return_code"]
            self.result.return_code = self._return_code


def attach_local(executor, run, *, callback=None, output_listener=None, start=True):
    if run.execution.endpoint.executor != "local":
        raise ValueError("Run does not belong to the local backend")
    status = detached_local_status(run)
    if status is None:
        raise NotImplementedError("PID-only local runs cannot be reattached")
    files = run.files
    spec = ExecutionSpec(
        argv=(),
        input_file=files.input_file,
        output_file=files.output_file,
        error_file=None if files.stderr is None else files.stderr.local_path,
        directory=run.execution.directory.local,
    )
    runner = LocalPersistentRunner(
        executor, spec, callback=callback, output_listener=output_listener,
        start=False, backend_id=run.execution.backend_id, state=status.state,
        failure_phase=status.failure_phase, failure_message=status.message,
    )
    runner._return_code = status.return_code
    runner.result.return_code = status.return_code
    if start:
        runner.start()
    return runner
