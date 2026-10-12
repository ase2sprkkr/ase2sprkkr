"""Local subprocess execution backend."""

from __future__ import annotations

import asyncio
import errno
import os
import sys
from dataclasses import replace
from threading import Event, Lock, Thread

from .base import ExecutionSpec, Executor, RunReference, RunState


def _mpi_command(launch):
    """Materialize portable MPI intent for this machine."""
    from ..configuration import mpi_runner

    command = mpi_runner(launch.mpi)
    if not command:
        if launch.mpi == "auto":
            return None
        raise FileNotFoundError(
            errno.ENOENT,
            "MPI was requested but no MPI launcher was found",
        )

    command = list(command)
    if launch.process_count is not None:
        launcher = os.path.basename(os.fspath(command[0]))
        count_option = "-n" if launcher == "srun" else "-np"
        command.extend((count_option, str(launch.process_count)))
    return tuple(command)


def materialize_local_spec(spec: ExecutionSpec) -> ExecutionSpec:
    """Choose serial/MPI argv immediately before local process creation."""
    if not spec.launch or spec.mpi_argv is None:
        return spec

    launcher = _mpi_command(spec.launch)
    if launcher is None:
        # ``auto`` may fall back to the serial command on this machine.
        return spec

    return replace(
        spec,
        argv=(*launcher, *spec.mpi_argv),
        stdin=None,
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

    The runner owns process lifecycle, cancellation, output-file writing and
    delivery of raw stdout/stderr to listeners. The output reader only parses
    streams and builds the TaskResult.
    """

    def __init__(self, reader, spec: ExecutionSpec, callback=None):
        self.reader = reader
        self.spec = spec
        self.callback = callback
        self.result = reader.create_result(spec.output_file, spec.input_file)

        self._state = RunState.NEW
        self._listeners = []
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

    @property
    def state(self):
        return self._state

    @property
    def reference(self):
        if self._pid is None:
            return None
        return RunReference(
            executor="local",
            native_id=str(self._pid),
            directory=os.path.abspath(self.spec.directory),
            input_file=_path(self.spec.input_file, self.spec.directory),
            output_file=_path(self.spec.output_file, self.spec.directory),
            error_file=_path(self.spec.error_file, self.spec.directory),
        )

    def add_output_listener(self, callback):
        """Receive raw output as ``callback(text, kind)``.

        ``kind`` is ``"out"`` or ``"err"``.
        """
        with self._lock:
            self._listeners.append(callback)

    def _emit_output(self, data, kind):
        if not data:
            return

        text = data.decode("utf8")
        if self.reader.print_output:
            print(text, end="", file=sys.stderr if kind == "err" else sys.stdout)

        with self._lock:
            listeners = tuple(self._listeners)
        for callback in listeners:
            callback(text, kind)

    def _output_target(self, kind):
        if kind == "err" and self.spec.error_file is not None:
            return self.spec.error_file
        return self.spec.output_file

    def _write_output(self, data, kind):
        target = self._output_target(kind)
        if target is not None:
            target.write(data)
            target.flush()

    async def _pump(self, source, target, kind):
        try:
            while data := await source.read(65536):
                self._write_output(data, kind)
                self._emit_output(data, kind)
                target.feed_data(data)
        finally:
            target.feed_eof()

    async def _stop_process(self):
        process = self._process
        if process is None or process.returncode is not None:
            return

        try:
            process.terminate()
        except ProcessLookupError:
            return

        try:
            await asyncio.wait_for(process.wait(), 0.1)
        except asyncio.TimeoutError:
            if process.returncode is not None:
                return
            try:
                process.kill()
            except ProcessLookupError:
                return
            await process.wait()

    def _ensure_stop_task(self):
        if self._process is None or self._process.returncode is not None:
            return None
        if self._stop_task is None or self._stop_task.done():
            self._stop_task = asyncio.create_task(self._stop_process())
        return self._stop_task

    def _schedule_stop(self):
        self._ensure_stop_task()

    def cancel(self):
        """Request cancellation from any thread."""
        self._cancel_requested.set()

        loop = self._loop
        if loop is None:
            if self._state == RunState.NEW:
                self._state = RunState.CANCELLED
            return

        try:
            loop.call_soon_threadsafe(self._schedule_stop)
        except RuntimeError:
            # The owning event loop has already stopped.
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
            raise RuntimeError("Run was cancelled before process start")

        try:
            self._process = await asyncio.create_subprocess_exec(
                *self.spec.argv,
                stdin=self.spec.stdin,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=self.spec.directory,
            )
            self._pid = self._process.pid
            self._state = RunState.RUNNING
            self._started.set()

            if self._cancel_requested.is_set():
                stop_task = self._ensure_stop_task()
                if stop_task is not None:
                    await stop_task

            stdout = asyncio.StreamReader()
            stderr = asyncio.StreamReader()

            tasks = [
                asyncio.create_task(
                    self.reader.parse_streams(stdout, stderr, self.result)
                ),
                asyncio.create_task(self._process.wait()),
                asyncio.create_task(self._pump(self._process.stdout, stdout, "out")),
                asyncio.create_task(self._pump(self._process.stderr, stderr, "err")),
            ]

            try:
                (output, error), return_code, _, _ = await asyncio.gather(*tasks)
            finally:
                for task in tasks:
                    if not task.done():
                        task.cancel()
                await asyncio.gather(*tasks, return_exceptions=True)

            result = self.reader.complete_result(
                self.result, output, error, return_code
            )
            self._state = RunState.FINISHED
            return result
        except BaseException:
            self._started.set()
            self._state = (
                RunState.CANCELLED
                if self._cancel_requested.is_set()
                else RunState.FAILED
            )
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
        finally:
            try:
                if self.callback:
                    self.callback(self.result)
            except BaseException as exc:
                if self._exception is None:
                    self._exception = exc
            self._done.set()

    def submit(self):
        """Start the run and return once its persistent process id is known."""
        with self._lock:
            if self._state == RunState.CANCELLED and self._thread is None:
                raise RuntimeError("Run was cancelled before submission")
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")

            self._state = RunState.QUEUED
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

        self._started.wait()

    def wait(self):
        """Wait for a submitted run and return its TaskResult."""
        if self._state == RunState.NEW:
            raise RuntimeError("Run has not been submitted")
        if self._state == RunState.CANCELLED and self._thread is None:
            raise RuntimeError("Run was cancelled before submission")

        self._done.wait()
        if self._exception is not None:
            raise self._exception
        return self.result

    def run(self):
        """Submit, if needed, and wait for completion."""
        if self._state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        """Async wrapper around submit/wait."""
        if self._state == RunState.NEW:
            self.submit()
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.wait)


class LocalExecutor(Executor):
    """Create runners executed as local subprocesses."""

    name = "local"

    def create_run(self, reader, spec: ExecutionSpec, callback=None):
        return LocalRunner(reader, materialize_local_spec(spec), callback=callback)
