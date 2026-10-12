"""PSI/J execution backend."""

from __future__ import annotations

import asyncio
import os
import sys
from pathlib import Path
from threading import Event, Lock, Thread

from .base import ExecutionSpec, Executor, RunReference, RunState


_STATE_MAP = {
    "NEW": RunState.NEW,
    "STAGE_IN": RunState.QUEUED,
    "QUEUED": RunState.QUEUED,
    "ACTIVE": RunState.RUNNING,
    "STAGE_OUT": RunState.RUNNING,
    "COMPLETED": RunState.FINISHED,
    "FAILED": RunState.FAILED,
    "CANCELED": RunState.CANCELLED,
}


def _load_psij():
    try:
        import psij
    except ImportError as exc:
        raise ImportError(
            "PSI/J support requires the optional 'psij-python' package. "
            "Install it with 'python -m pip install \"ase2sprkkr[psij]\"'."
        ) from exc
    return psij


def _path(value, directory=None):
    if value is None:
        return None
    value = getattr(value, "name", value)
    try:
        path = os.fspath(value)
    except TypeError as exc:
        raise ValueError("PSI/J execution requires named input/output files") from exc
    if directory and not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


class PsijRunner:
    """One job managed by a PSI/J executor.

    PSI/J scheduler executors redirect stdout/stderr to files. The runner
    follows those files while the job is active and feeds the same asynchronous
    output reader that is used by :class:`LocalRunner`.
    """

    def __init__(
        self,
        executor,
        reader,
        spec: ExecutionSpec,
        callback=None,
        *,
        job=None,
        attached=False,
    ):
        self.executor = executor
        self.reader = reader
        self.spec = spec
        self.callback = callback
        self.result = reader.create_result(spec.output_file, spec.input_file)

        self._stdout_path = _path(spec.output_file, spec.directory)
        if self._stdout_path is None:
            raise ValueError("PSI/J execution requires an output file")
        self._stderr_path = _path(spec.error_file, spec.directory) or self._stdout_path + ".err"
        self._stdin_path = _path(spec.stdin, spec.directory)

        self._listeners = []
        self._lock = Lock()
        self._done = Event()
        self._exception = None
        self._thread = None
        self._attached = attached

        self._psij = executor._psij
        self._job = job or self._create_job()

        if attached:
            self._start_monitor()

    def _mpi_launcher(self):
        # Slurm's native launcher knows the scheduler allocation. Other PSI/J
        # backends commonly use mpirun; users needing a custom command can use
        # an explicit MPI launcher in the profile instead.
        return "srun" if self.executor.backend.lower() == "slurm" else "mpirun"

    def _job_command(self):
        launch = self.spec.launch
        if not launch or self.spec.mpi_argv is None:
            return tuple(self.spec.argv), None, None

        if isinstance(launch.mpi, tuple):
            # An explicit launcher is already a concrete command. Preserve it
            # instead of asking PSI/J to reinterpret it.
            return (*launch.mpi, *self.spec.mpi_argv), None, None

        resources = None
        if launch.process_count is not None:
            resources = self._psij.ResourceSpecV1(
                process_count=launch.process_count
            )
        return tuple(self.spec.mpi_argv), resources, self._mpi_launcher()

    def _create_job(self):
        argv, resources, launcher = self._job_command()
        if not argv:
            raise ValueError("Cannot submit a PSI/J run without a command")
        job_spec = self._psij.JobSpec(
            executable=os.fspath(argv[0]),
            arguments=[os.fspath(arg) for arg in argv[1:]],
            directory=self.spec.directory,
            stdin_path=None if self.spec.launch and self.spec.mpi_argv is not None else self._stdin_path,
            stdout_path=self._stdout_path,
            stderr_path=self._stderr_path,
            resources=resources,
            launcher=launcher,
        )
        return self._psij.Job(job_spec)

    @property
    def state(self):
        state = getattr(getattr(self._job, "status", None), "state", None)
        return _STATE_MAP.get(getattr(state, "name", None), RunState.UNKNOWN)

    @property
    def reference(self):
        native_id = getattr(self._job, "native_id", None)
        if native_id is None:
            return None
        return RunReference(
            executor=f"psij:{self.executor.backend}",
            native_id=str(native_id),
            computer=self.executor.url,
            directory=os.fspath(self.spec.directory),
            input_file=_path(self.spec.input_file, self.spec.directory),
            output_file=self._stdout_path,
            error_file=self._stderr_path,
        )

    def add_output_listener(self, callback):
        """Receive raw output as ``callback(text, kind)``."""
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

    @staticmethod
    def _flush(value):
        if value is not None and hasattr(value, "flush"):
            value.flush()

    @staticmethod
    def _close(value):
        if value is not None and hasattr(value, "close"):
            value.close()

    def _prepare_io(self):
        self._flush(self.spec.input_file)
        self._flush(self.spec.output_file)
        self._flush(self.spec.error_file)

        # PSI/J opens stdout/stderr itself. Do not keep another writer open.
        self._close(self.spec.output_file)
        self._close(self.spec.error_file)

        Path(self._stdout_path).parent.mkdir(parents=True, exist_ok=True)
        Path(self._stderr_path).parent.mkdir(parents=True, exist_ok=True)
        Path(self._stdout_path).write_bytes(b"")
        Path(self._stderr_path).write_bytes(b"")

    async def _follow_file(self, path, stream, kind, finished):
        position = 0
        empty_after_finish = False
        try:
            while True:
                data = b""
                try:
                    size = os.path.getsize(path)
                    if size < position:
                        position = 0
                    with open(path, "rb") as source:
                        source.seek(position)
                        data = source.read()
                except FileNotFoundError:
                    pass

                if data:
                    position += len(data)
                    self._emit_output(data, kind)
                    stream.feed_data(data)
                    empty_after_finish = False
                elif finished.is_set():
                    if empty_after_finish:
                        break
                    empty_after_finish = True

                await asyncio.sleep(self.executor.poll_interval)
        finally:
            stream.feed_eof()

    async def _monitor(self):
        stdout = asyncio.StreamReader()
        stderr = asyncio.StreamReader()
        finished = asyncio.Event()

        parse_task = asyncio.create_task(
            self.reader.parse_streams(stdout, stderr, self.result)
        )
        followers = [
            asyncio.create_task(self._follow_file(self._stdout_path, stdout, "out", finished)),
            asyncio.create_task(self._follow_file(self._stderr_path, stderr, "err", finished)),
        ]

        status = None
        try:
            status = await asyncio.to_thread(self._job.wait)
        finally:
            finished.set()
            await asyncio.gather(*followers, return_exceptions=False)

        output, error = await parse_task
        return_code = getattr(status, "exit_code", None)
        if return_code is None:
            return_code = 0 if getattr(getattr(status, "state", None), "name", None) == "COMPLETED" else 1
        return self.reader.complete_result(self.result, output, error, return_code)

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
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

    def _start_monitor(self):
        with self._lock:
            if self._thread is not None:
                return
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def submit(self):
        """Submit the job through PSI/J and return after the scheduler accepts it."""
        if self._thread is not None or self.state != RunState.NEW:
            raise RuntimeError(f"Run is already {self.state.value}")
        self._prepare_io()
        self.executor._job_executor.submit(self._job)
        self._start_monitor()

    def wait(self):
        """Wait for the PSI/J job and output parsing to finish."""
        if self._thread is None:
            raise RuntimeError("Run has not been submitted")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        return self.result

    def run(self):
        if self._thread is None:
            self.submit()
        return self.wait()

    async def run_async(self):
        if self._thread is None:
            self.submit()
        return await asyncio.to_thread(self.wait)

    def cancel(self):
        """Cancel the scheduler job through PSI/J."""
        if self._thread is None and not self._attached:
            raise RuntimeError("Run has not been submitted")
        self._job.cancel()


class PsijExecutor(Executor):
    """Execute calculations through a PSI/J backend such as Slurm or PBS."""

    name = "psij"

    def __init__(
        self,
        backend="slurm",
        *,
        url=None,
        config=None,
        job_executor=None,
        poll_interval=0.1,
    ):
        self.backend = backend
        self.url = url
        self.poll_interval = poll_interval
        self._psij = _load_psij()
        self._job_executor = (
            job_executor
            if job_executor is not None
            else self._psij.JobExecutor.get_instance(backend, url=url, config=config)
        )

    def create_run(self, reader, spec: ExecutionSpec, callback=None):
        return PsijRunner(self, reader, spec, callback=callback)

    def attach(self, reference: RunReference, reader, callback=None):
        expected = f"psij:{self.backend}"
        if reference.executor != expected:
            raise ValueError(
                f"Run reference belongs to {reference.executor!r}, not {expected!r}"
            )
        if not reference.directory or not reference.output_file:
            raise ValueError("PSI/J run reference lacks output location")

        job = self._psij.Job()
        self._job_executor.attach(job, reference.native_id)
        spec = ExecutionSpec(
            argv=(),
            input_file=reference.input_file,
            output_file=reference.output_file,
            error_file=reference.error_file,
            directory=reference.directory,
        )
        return PsijRunner(
            self,
            reader,
            spec,
            callback=callback,
            job=job,
            attached=True,
        )
