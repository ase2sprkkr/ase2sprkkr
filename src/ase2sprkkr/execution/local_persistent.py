"""Reconnectable local execution using the detached execution helper."""

from __future__ import annotations

import asyncio
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path
from threading import Event, Lock, Thread

from .base import (
    ExecutionSpec,
    FailurePhase,
    OutputObservers,
    RunReference,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .detached import (
    _request_path,
    _run_directory,
    cancel_detached,
    detached_failure_phase,
    detached_launch_failed,
    detached_run_state,
    detached_stream_file,
    status_detached,
    submit_detached,
)
from .outcome import ExecutionOutcome
from .transport import AtomicLocalFile

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


def _reference(run_id, spec):
    directory = os.path.abspath(os.fspath(spec.directory))
    return RunReference(
        executor="local",
        native_id=run_id,
        directory=directory,
        input_file=_file_path(spec.input_file, directory),
        output_file=_file_path(spec.output_file, directory),
        error_file=_file_path(spec.error_file, directory),
    )


def _reference_from_request(run_id, directory, request):
    return RunReference(
        executor="local",
        native_id=run_id,
        directory=os.path.abspath(os.fspath(directory)),
        input_file=request.get("client_input_file") or request.get("input_file"),
        output_file=request.get("client_output_file"),
        error_file=request.get("client_error_file"),
    )


def detached_local_status(reference):
    """Return detached-helper status, or ``None`` for legacy PID references."""
    if not reference.directory:
        return None
    try:
        state = status_detached(reference.directory, reference.native_id)
    except (FileNotFoundError, json.JSONDecodeError):
        return None
    return RunStatus(
        reference,
        detached_run_state(state.get("state", "unknown")),
        return_code=state.get("return_code"),
        message=state.get("error"),
        failure_phase=detached_failure_phase(state),
    )


def recover_local_submission(run_id, submission):
    """Recover a local detached submit whose reference was not persisted."""
    directory = os.path.abspath(os.fspath(submission.directory))
    try:
        state = status_detached(directory, run_id)
        request = json.loads(_request_path(directory, run_id).read_text())
    except FileNotFoundError:
        return None
    reference = _reference_from_request(run_id, directory, request)
    return RunStatus(
        reference,
        detached_run_state(state.get("state", "unknown")),
        return_code=state.get("return_code"),
        message=state.get("error"),
        failure_phase=detached_failure_phase(state),
    )


class LocalPersistentRunner:
    """Detached local execution without task-specific output interpretation."""

    def __init__(
        self,
        executor,
        spec: ExecutionSpec,
        *,
        run_id=None,
        callback=None,
        output_listener=None,
        print_output=False,
        start=False,
        reference=None,
        state=RunState.NEW,
        failure_phase=None,
        failure_message=None,
    ):
        self.executor = executor
        self.spec = spec
        self.callback = callback
        self.print_output = bool(print_output)
        self.result = ExecutionOutcome(spec.output_file, spec.input_file)

        self._run_id = run_id or (reference.native_id if reference is not None else None)
        self._reference = reference
        self._attached = reference is not None
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
    def reference(self):
        return self._reference

    @property
    def state(self):
        if self._reference is None or self._state in _FINAL_STATES:
            return self._state
        status = detached_local_status(self._reference)
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
        if not data:
            return
        text = data.decode("utf8")
        if self.print_output:
            print(text, end="", file=sys.stderr if kind == "err" else sys.stdout)
        self._output_observers.emit(text, kind)

    def _open_target(self, path):
        if path is None:
            return None
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if self._attached:
            return AtomicLocalFile(path)
        return path.open("wb")

    def _open_mirrors(self):
        output_path = _file_path(self.spec.output_file, self.spec.directory)
        error_path = _file_path(self.spec.error_file, self.spec.directory)
        output = self._open_target(output_path)
        if error_path is None or error_path == output_path:
            error = output
        else:
            error = self._open_target(error_path)
        return output, error

    @staticmethod
    def _finish_target(target, *, commit):
        if target is None:
            return
        method = getattr(target, "commit" if commit else "discard", None)
        if method is not None:
            method()
        else:
            target.close()

    @classmethod
    def _finish_mirrors(cls, output, error, *, commit):
        if error is not None and error is not output:
            cls._finish_target(error, commit=commit)
        cls._finish_target(output, commit=commit)

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
        output_target = error_target = None
        mirrors_committed = False
        started_streams = False
        try:
            while True:
                status = status_detached(
                    self._reference.directory, self._reference.native_id
                )
                self._failure_phase = detached_failure_phase(status)
                self._failure_message = status.get("error")
                if detached_launch_failed(status):
                    self._state = RunState.FAILED
                    self._return_code = None
                    return self.result

                state = detached_run_state(status.get("state", "unknown"))
                started = status.get("started_at") is not None
                if started and not started_streams:
                    output_target, error_target = self._open_mirrors()
                    started_streams = True

                had_data = False
                if started:
                    for kind, key, target in (
                        ("out", "output_file", output_target),
                        ("err", "error_file", error_target),
                    ):
                        data = self._read_new(
                            detached_stream_file(status, key), offsets[kind]
                        )
                        if not data:
                            continue
                        offsets[kind] += len(data)
                        if target is not None:
                            target.write(data)
                            target.flush()
                        self._emit_output(data, kind)
                        had_data = True

                if state in _FINAL_STATES and not had_data:
                    self._state = state
                    self._return_code = status.get("return_code")
                    if self._return_code is None and started:
                        self._return_code = 0 if state == RunState.FINISHED else 1
                    self.result.return_code = self._return_code
                    if started_streams:
                        self._finish_mirrors(
                            output_target, error_target, commit=True
                        )
                        mirrors_committed = True
                    return self.result
                if not had_data:
                    await asyncio.sleep(0.05)
        finally:
            if started_streams and not mirrors_committed:
                self._finish_mirrors(output_target, error_target, commit=False)

    def _close_output_files(self):
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
            self._close_output_files()
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
            if self._reference is None:
                raise RuntimeError("Run has no backend reference")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def start(self):
        if self._reference is None:
            return self.submit()
        return self._start_monitor()

    def _request(self):
        directory = os.path.abspath(os.fspath(self.spec.directory))
        run_directory = _run_directory(directory, self._run_id)
        return {
            "id": self._run_id,
            "executor": "local",
            "directory": directory,
            "argv": list(self.spec.argv),
            "mpi_argv": None
            if self.spec.mpi_argv is None
            else list(self.spec.mpi_argv),
            "launch": asdict(self.spec.launch),
            "input_file": _file_path(self.spec.input_file, directory),
            "stdin": _file_path(self.spec.stdin, directory),
            "output_file": os.fspath(run_directory / "stdout"),
            "error_file": os.fspath(run_directory / "stderr"),
            "client_input_file": _file_path(self.spec.input_file, directory),
            "client_output_file": _file_path(self.spec.output_file, directory),
            "client_error_file": _file_path(self.spec.error_file, directory),
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
            self._reference = _reference(self._run_id, self.spec)
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
        if self._reference is None:
            raise RuntimeError("Run has not been submitted")
        state = cancel_detached(
            self._reference.directory, self._reference.native_id
        )
        self._state = detached_run_state(state.get("state", "cancelled"))
        self._failure_phase = detached_failure_phase(state)
        self._failure_message = state.get("error")
        if state.get("return_code") is not None:
            self._return_code = state["return_code"]
            self.result.return_code = self._return_code


def attach_local(
    executor,
    reference,
    *,
    callback=None,
    output_listener=None,
    start=True,
    print_output=False,
):
    if reference.executor != "local":
        raise ValueError(
            f"Run reference belongs to {reference.executor!r}, not 'local'"
        )
    status = detached_local_status(reference)
    if status is None:
        raise NotImplementedError("Legacy PID-only local runs cannot be reattached")
    spec = ExecutionSpec(
        argv=(),
        input_file=reference.input_file,
        output_file=reference.output_file,
        error_file=reference.error_file,
        directory=reference.directory or ".",
    )
    runner = LocalPersistentRunner(
        executor,
        spec,
        callback=callback,
        output_listener=output_listener,
        print_output=print_output,
        start=False,
        reference=reference,
        state=status.state,
        failure_phase=status.failure_phase,
        failure_message=status.message,
    )
    runner._return_code = status.return_code
    runner.result.return_code = status.return_code
    if start:
        runner.start()
    return runner
