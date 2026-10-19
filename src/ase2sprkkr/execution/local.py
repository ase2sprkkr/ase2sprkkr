"""Local subprocess execution backend."""

from __future__ import annotations

import asyncio
import os
from threading import Event, Lock, Thread

from .base import (
    ExecutionOutcome,
    ExecutionSpec,
    Executor,
    FailurePhase,
    OutputObservers,
    finish_monitor,
    wait_for_monitor,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .file_availability import clear_declared_outputs
from .local_process import (
    materialize_local_spec,
    pump_stream,
    start_local_process,
    stop_local_process,
)


def _path(value, directory):
    if value is None:
        return None
    value = getattr(value, "name", value)
    if isinstance(value, int):
        return None
    try:
        path = os.fspath(value)
    except TypeError:
        return None
    if not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


class LocalRunner:
    """One local subprocess execution.

    The runner owns only process lifecycle, output-file writing and raw output
    observation. It deliberately performs no task-specific parsing.
    """

    def __init__(self, spec: ExecutionSpec, callback=None):
        self.spec = spec
        self.callback = callback
        self.result = ExecutionOutcome(
            output_file=_path(spec.output_file, spec.directory),
            input_file=_path(spec.input_file, spec.directory),
            error_file=_path(spec.error_file, spec.directory),
        )

        self._state = RunState.NEW
        self._output_observers = OutputObservers()
        self._lock = Lock()
        self._cancel_requested = Event()
        self._started = Event()
        self._done = Event()

        self._thread = None
        self._loop = None
        self._process = None
        self._pid = None
        self._stop_task = None
        self._exception = None
        self._callback_exception = None
        self._failure_phase = None
        self._failure_message = None

    @property
    def state(self):
        return self._state

    @property
    def exception(self):
        return self._exception

    @property
    def backend_id(self):
        return None if self._pid is None else str(self._pid)

    @property
    def failure_phase(self):
        return self._failure_phase

    @property
    def failure_message(self):
        return self._failure_message

    @property
    def return_code(self):
        return self.result.return_code

    def add_output_listener(self, callback):
        self._output_observers.add(callback)

    def add_output_error_listener(self, callback):
        self._output_observers.add_error_listener(callback)

    def _emit_output(self, data, kind):
        if data:
            self._output_observers.emit(data.decode("utf8", errors="replace"), kind)

    def _output_target(self, kind):
        # A missing stderr target means live observation only, not an
        # implicit merge into the saved stdout artifact.
        return self.spec.error_file if kind == "err" else self.spec.output_file

    def _forward_output(self, data, kind):
        target = self._output_target(kind)
        if target is not None:
            target.write(data)
            target.flush()
        self._emit_output(data, kind)

    async def _stop_process(self):
        await stop_local_process(self._process)

    def _ensure_stop_task(self):
        if self._process is None or self._process.returncode is not None:
            return None
        if self._stop_task is None or self._stop_task.done():
            self._stop_task = asyncio.create_task(self._stop_process())
        return self._stop_task

    def _schedule_stop(self):
        self._ensure_stop_task()

    def cancel(self):
        self._cancel_requested.set()
        loop = self._loop
        if loop is None:
            if self._state == RunState.NEW:
                self._state = RunState.CANCELLED
            return
        try:
            loop.call_soon_threadsafe(self._schedule_stop)
        except RuntimeError:
            pass

    def _close_output_files(self):
        closed = set()
        for file in (self.spec.output_file, self.spec.error_file):
            if file is None or id(file) in closed:
                continue
            closed.add(id(file))
            if hasattr(file, "close"):
                file.close()

    async def _execute(self):
        self._loop = asyncio.get_running_loop()
        if self._cancel_requested.is_set():
            self._state = RunState.CANCELLED
            self._started.set()
            return self.result

        try:
            protected = [
                _path(value, self.spec.directory)
                for value in (self.spec.input_file, self.spec.stdin,
                              *self.spec.input_files, *self.spec.cache_files)
            ]
            # Preserve already-open stdout/stderr handles: an unlink would
            # disconnect their inode from the named destination.
            outputs = [
                _path(value, self.spec.directory)
                for value in self.spec.output_files
            ]
            clear_declared_outputs(
                outputs,
                protected=(*protected, self.result.output_file, self.result.error_file),
            )
            try:
                self._process = await start_local_process(self.spec)
            except BaseException as exc:
                self._failure_phase = FailurePhase.LAUNCH
                self._failure_message = str(exc)
                self._state = RunState.FAILED
                raise

            self._pid = self._process.pid
            self._state = RunState.RUNNING
            self._started.set()

            if self._cancel_requested.is_set():
                stop_task = self._ensure_stop_task()
                if stop_task is not None:
                    await stop_task

            wait_task = asyncio.create_task(self._process.wait())
            stdout_task = asyncio.create_task(
                pump_stream(
                    self._process.stdout,
                    lambda data: self._forward_output(data, "out"),
                )
            )
            stderr_task = asyncio.create_task(
                pump_stream(
                    self._process.stderr,
                    lambda data: self._forward_output(data, "err"),
                )
            )
            tasks = (wait_task, stdout_task, stderr_task)
            try:
                return_code, _, _ = await asyncio.gather(*tasks)
            finally:
                for task in tasks:
                    if not task.done():
                        task.cancel()
                await asyncio.gather(*tasks, return_exceptions=True)

            self.result.return_code = return_code
            if self._cancel_requested.is_set():
                self._state = RunState.CANCELLED
            elif return_code == 0:
                self._state = RunState.FINISHED
            else:
                self._state = RunState.FAILED
                self._failure_phase = FailurePhase.EXECUTION
                self._failure_message = f"Process exited with status {return_code}"
            return self.result
        except BaseException:
            self._started.set()
            if self._state not in (RunState.CANCELLED, RunState.FAILED):
                self._state = RunState.FAILED
                self._failure_phase = FailurePhase.EXECUTION
            if self._process is not None and self._process.returncode is None:
                await self._stop_process()
            raise
        finally:
            if self._stop_task is not None:
                await asyncio.gather(self._stop_task, return_exceptions=True)
            self._close_output_files()
            self._process = None
            self._loop = None

    def _thread_main(self):
        try:
            asyncio.run(self._execute())
        except BaseException as exc:
            self._exception = exc
            if self._failure_message is None:
                self._failure_message = str(exc)
        finally:
            finish_monitor(self)

    def submit(self):
        with self._lock:
            if self._state == RunState.CANCELLED and self._thread is None:
                raise RuntimeError("Run was cancelled before submission")
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")
            self._state = RunState.SUBMITTING
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

        self._started.wait()
        if self._pid is None and self._state == RunState.FAILED:
            self._done.wait()
            error = self._exception or RuntimeError(self._failure_message or "Launch failed")
            raise SubmissionNotStarted(str(error)) from self._exception

    def wait(self):
        if self._state == RunState.NEW:
            raise RuntimeError("Run has not been submitted")
        if self._state == RunState.CANCELLED and self._thread is None:
            return self.result
        return wait_for_monitor(self)

    def run(self):
        if self._state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        if self._state == RunState.NEW:
            self.submit()
        return await asyncio.to_thread(self.wait)


class LocalExecutor(Executor):
    """Create raw runners executed as local subprocesses."""

    name = "local"

    def create_run(self, spec: ExecutionSpec, callback=None):
        return LocalRunner(materialize_local_spec(spec), callback=callback)

    def create_project_run(
        self, spec: ExecutionSpec, *, run_id: str, callback=None
    ):
        from .local_project import LocalProjectRunner

        return LocalProjectRunner(
            self,
            spec,
            run_id=run_id,
            callback=callback,
        )

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        from .remote_service import recover_remote
        from .local_project import local_manifest
        reply = recover_remote(run_id, {
            "executor": "local",
            "directory": submission.execution.directory.local,
        })
        if reply is None:
            # The absence of a ticket is not proof that submission never
            # happened. Never authorize an automatic duplicate launch.
            raise NotImplementedError("Submission identity is not recoverable")
        status = RunStatus(
            state=RunState(reply["state"]),
            backend_id=reply["remote_id"],
            return_code=reply.get("return_code"),
            message=reply.get("error"),
            failure_phase=reply.get("failure_phase"),
        )
        if status.state in (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED):
            from dataclasses import replace
            status = replace(status, files=local_manifest(submission.files, status))
        return status

    def attach(self, run, callback=None, output_listener=None, start=True):
        from .local_project import attach_local
        return attach_local(
            self, run, callback=callback,
            output_listener=output_listener, start=start,
        )

    def acknowledge_submission(self, run):
        from .remote_service import acknowledge_submission
        acknowledge_submission(run.id)

    def acknowledge(self, run):
        from .remote_service import acknowledge_remote
        acknowledge_remote(run.execution.backend_id)

    def status(self, *runs):
        from .local_project import local_status
        statuses = []
        for run in runs:
            if run.execution.endpoint.executor != "local":
                raise ValueError("Run does not belong to the local backend")
            statuses.append(local_status(run.execution.backend_id, run.files))
        return tuple(statuses)
