"""Raw reconnectable local process for a project-associated calculation."""

from __future__ import annotations

import asyncio
import os
from dataclasses import asdict
from threading import Event, Lock, Thread

from . import local_backend
from .growing_file_tail import GrowingFileTail
from .base import (
    ExecutionOutcome, ExecutionSpec, FailurePhase, OutputObservers,
    RunState, SubmissionNotStarted,
)


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


class LocalProjectRunner:
    """Client-side monitor of a worker identified by PID and birth timestamp."""

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
        if self._backend_id is None:
            return self._state
        status = local_backend.query(self.spec.directory, self._backend_id)
        self._state = status.state
        self._failure_phase = status.failure_phase
        self._failure_message = status.message
        self._return_code = status.return_code
        self.result.return_code = self._return_code
        if status.state in (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED):
            if self._files is not None:
                self._files = local_backend.local_manifest(
                    self._files, self.spec.directory, self._backend_id,
                )
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

    async def _monitor(self):
        offsets = {"out": 0, "err": 0}
        initial = local_backend.observed_paths(
            self.spec.directory, self._backend_id,
            self.result.output_file, self.result.error_file,
        )
        async with GrowingFileTail(
            dict(zip(("out", "err"), initial)), offsets=offsets,
            watch_paths=tuple(path for path in (
                self.result.output_file, self.result.error_file,
            ) if path is not None),
        ) as tail:
            while True:
                state = self._refresh()
                paths = local_backend.observed_paths(
                    self.spec.directory, self._backend_id,
                    self.result.output_file, self.result.error_file,
                )
                chunks = tail.read_available(dict(zip(("out", "err"), paths)))
                for kind, data in chunks:
                    self._output_observers.emit(
                        data.decode("utf8", errors="replace"), kind,
                    )
                if state in (RunState.FINISHED, RunState.FAILED,
                             RunState.CANCELLED, RunState.UNKNOWN) and not chunks:
                    return self.result
                if not chunks:
                    await tail.wait_for_change()

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
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
            self._failure_message = str(exc)
            self._failure_phase = self._failure_phase or FailurePhase.EXECUTION
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
                raise RuntimeError("Local run has no backend identity")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def _request(self):
        return {
            "argv": list(self.spec.argv),
            "input_file": self.result.input_file,
            "mpi_argv": None if self.spec.mpi_argv is None else list(self.spec.mpi_argv),
            "launch": asdict(self.spec.launch),
            "stdin": _file_path(self.spec.stdin, self.spec.directory),
            "output_file": self.result.output_file,
            "error_file": self.result.error_file,
        }

    def submit(self):
        if self._state != RunState.NEW:
            raise RuntimeError(f"Run is already {self._state.value}")
        if not self._run_id:
            raise RuntimeError("Project run has no run id")
        try:
            self._backend_id = local_backend.submit(
                self.spec.directory, self._run_id, self._request(),
            )
        except SubmissionNotStarted as exc:
            self._state = RunState.FAILED
            self._failure_phase = FailurePhase.SUBMISSION
            self._failure_message = str(exc)
            raise
        self._state = RunState.RUNNING
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
        local_backend.cancel(self.spec.directory, self._backend_id)


def attach_local(executor, run, *, callback=None, output_listener=None, start=True):
    if run.execution.endpoint.executor != "local":
        raise ValueError("Run does not belong to the local backend")
    backend_id = run.execution.backend_id
    if local_backend._parse_id(backend_id) is None:
        raise NotImplementedError("Cannot reattach legacy PID-only local runs")
    files = run.files
    spec = ExecutionSpec(
        argv=(), input_file=files.input_file,
        output_file=files.output_file,
        error_file=None if files.stderr is None else files.stderr.local_path,
        directory=run.execution.directory.local,
    )
    status = local_backend.query(spec.directory, backend_id)
    runner = LocalProjectRunner(
        executor, spec, callback=callback, output_listener=output_listener,
        backend_id=backend_id, run_id=run.id,
        state=status.state, failure_phase=status.failure_phase,
        failure_message=status.message,
    )
    runner.set_manifest(files)
    runner._return_code = status.return_code
    runner.result.return_code = status.return_code
    if start:
        runner.start()
    return runner
