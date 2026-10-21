"""PSI/J execution backend."""

from __future__ import annotations

import asyncio
import logging
from collections import OrderedDict
import os
import shlex
import subprocess
import time
from datetime import datetime, timedelta
from dataclasses import replace
from pathlib import Path
from threading import Event, Lock, RLock

from .file_availability import discover_local_files, clear_declared_outputs, validate_output_collisions
from .growing_file_tail import follow_growing_files
from .base import (
    ExecutionOutcome,
    ExecutionSpec,
    Executor,
    FailurePhase,
    OutputObservers,
    close_output_files,
    file_path,
    finish_monitor,
    wait_for_monitor,
    RunState,
    RunStatus,
    SubmissionNotStarted,
    request_cancellation,
    start_monitor,
)


_log = logging.getLogger(__name__)

_OBSERVATION_STOPPED = object()


def _ignore_job_status(job, status):
    """PSI/J 0.9.11 wraps None as a broken callable on unsubscribe."""
    pass


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

_SLURM_STATE_MAP = {
    "PENDING": RunState.QUEUED,
    "CONFIGURING": RunState.QUEUED,
    "SUSPENDED": RunState.QUEUED,
    "RUNNING": RunState.RUNNING,
    "COMPLETING": RunState.RUNNING,
    "RESIZING": RunState.RUNNING,
    "COMPLETED": RunState.FINISHED,
    "CANCELLED": RunState.CANCELLED,
    "CANCELED": RunState.CANCELLED,
    "BOOT_FAIL": RunState.FAILED,
    "DEADLINE": RunState.FAILED,
    "FAILED": RunState.FAILED,
    "NODE_FAIL": RunState.FAILED,
    "OUT_OF_MEMORY": RunState.FAILED,
    "PREEMPTED": RunState.FAILED,
    "REVOKED": RunState.FAILED,
    "TIMEOUT": RunState.FAILED,
}


def _psij_state(status, *, attached=False, return_code=None):
    """Interpret native PSI/J job facts consistently across all observers.

    A newly attached job reports NEW before PSI/J has received its real
    scheduler status. That is UNKNOWN, not evidence of a new submission.
    A completed job with a nonzero exit code is an execution failure.
    """
    name = getattr(getattr(status, "state", None), "name", None)
    if attached and name == "NEW":
        return RunState.UNKNOWN
    state = _STATE_MAP.get(name, RunState.UNKNOWN)
    if return_code is None:
        return_code = getattr(status, "exit_code", None)
    if state == RunState.FINISHED and return_code not in (None, 0):
        return RunState.FAILED
    return state


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
    try:
        return file_path(value, directory, required=True)
    except ValueError as exc:
        raise ValueError("PSI/J execution requires named input/output files") from exc


def _scheduler_output(argv):
    return subprocess.run(
        tuple(argv),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout


def _slurm_state(value):
    value = (value or "").strip().upper().split("+", 1)[0]
    return _SLURM_STATE_MAP.get(value, RunState.UNKNOWN)


def _slurm_exit_code(value):
    if not value:
        return None
    code = value.split(":", 1)[0]
    try:
        return int(code)
    except ValueError:
        return None


def _slurm_details(native_id):
    output = _scheduler_output(("scontrol", "--oneliner", "show", "job", str(native_id)))
    fields = {}
    for token in shlex.split(output):
        name, separator, value = token.partition("=")
        if separator:
            fields[name] = value
    return fields


def _slurm_current_matches(run_id):
    marker = f"ase2sprkkr:{run_id}"
    output = _scheduler_output(
        ("squeue", "--noheader", "--me", "--format=%i|%k|%T|%Z")
    )
    matches = []
    for line in output.splitlines():
        fields = line.split("|", 3)
        if len(fields) != 4:
            continue
        native_id, comment, state, directory = (field.strip() for field in fields)
        if comment == marker:
            matches.append((native_id, state, directory))
    return tuple(matches)


def _slurm_accounting_stores_comments():
    output = _scheduler_output(("scontrol", "show", "config"))
    for line in output.splitlines():
        name, separator, value = line.partition("=")
        if separator and name.strip() == "AccountingStoreFlags":
            flags = {flag.strip().lower() for flag in value.split(",")}
            return "job_comment" in flags
    return False


def _slurm_time(value):
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is not None:
        parsed = parsed.astimezone()
    return parsed.strftime("%Y-%m-%dT%H:%M:%S")


def _slurm_accounting_matches(run_id, submitted_at):
    marker = f"ase2sprkkr:{run_id}"
    output = _scheduler_output(
        (
            "sacct",
            "-X",
            "-n",
            "-P",
            f"--starttime={_slurm_time(submitted_at)}",
            "--format=JobIDRaw,Comment,State,WorkDir,StdOut,StdErr,ExitCode",
        )
    )
    matches = []
    for line in output.splitlines():
        fields = line.rstrip("|").split("|")
        if len(fields) != 7:
            continue
        native_id, comment, state, directory, stdout, stderr, exit_code = fields
        if comment.strip() == marker:
            matches.append(
                (
                    native_id.strip(),
                    state.strip(),
                    directory.strip(),
                    stdout.strip() or None,
                    stderr.strip() or None,
                    _slurm_exit_code(exit_code.strip()),
                )
            )
    return tuple(matches)


def _slurm_status_native(native_id):
    """Read scheduler facts for a local Slurm job without initializing PSI/J."""
    native_id = str(native_id)
    # A PSI/J attach initially reports NEW; short-lived remote commands
    # cannot wait for an asynchronous update or retain the Job cache.
    # Query the same scheduler facts for local and SSH PSI/J runs.
    try:
        details = _slurm_details(native_id)
    except (FileNotFoundError, subprocess.CalledProcessError):
        details = {}

    if details:
        state = _slurm_state(details.get("JobState"))
        return_code = _slurm_exit_code(details.get("ExitCode"))
    else:
        try:
            output = _scheduler_output((
                "sacct", "-X", "-n", "-P", "-j", native_id,
                "--format=JobIDRaw,State,ExitCode",
            ))
        except (FileNotFoundError, subprocess.CalledProcessError):
            output = ""
        state = RunState.UNKNOWN
        return_code = None
        for line in output.splitlines():
            fields = line.rstrip("|").split("|")
            if len(fields) == 3 and fields[0].strip() == native_id:
                state = _slurm_state(fields[1])
                return_code = _slurm_exit_code(fields[2].strip())
                break

    if not state.is_terminal:
        return_code = None
    elif state == RunState.FINISHED and return_code is None:
        return_code = 0
    elif state == RunState.FINISHED and return_code != 0:
        state = RunState.FAILED
    return RunStatus(
        state=state,
        backend_id=native_id,
        return_code=return_code,
        failure_phase=FailurePhase.EXECUTION if state == RunState.FAILED else None,
    )



class PsijRunner:
    """One raw job managed by a PSI/J executor."""

    def __init__(
        self,
        executor,
        spec: ExecutionSpec,
        callback=None,
        *,
        output_listener=None,
        start=False,
        job=None,
        attached=False,
        run_id=None,
    ):
        self.executor = executor
        self.spec = spec
        self.callback = callback
        self._run_id = run_id

        self._stdout_path = _path(spec.output_file, spec.directory)
        self._stderr_path = _path(spec.error_file, spec.directory)
        self._stdin_path = _path(spec.stdin, spec.directory)
        self.result = ExecutionOutcome(
            output_file=self._stdout_path,
            input_file=_path(spec.input_file, spec.directory),
            error_file=self._stderr_path,
        )

        self._output_observers = OutputObservers(output_listener)
        self._lock = Lock()
        self._submit_lock = Lock()
        self._done = Event()
        self._exception = None
        self._callback_exception = None
        self._return_code = None
        self._failure_phase = None
        self._failure_message = None
        self._thread = None
        self._attached = attached
        self._files = None
        self._cancel_future = None
        self._observed_state = None
        self._observation_stopped = Event()
        self._observation_signal_lock = Lock()
        self._observation_wakeup = None

        self._psij = executor._psij
        self._job = job or self._create_job()

        if start:
            self.start()

    def _mpi_launcher(self):
        return "srun" if self.executor.backend.lower() == "slurm" else "mpirun"

    def _job_command(self):
        launch = self.spec.launch
        if not launch or self.spec.mpi_argv is None:
            return tuple(self.spec.argv), None, None
        mpi = launch.mpi
        if mpi.launcher:
            return (*mpi.launcher, *self.spec.mpi_argv), None, None
        resources = None
        if mpi.process_count is not None:
            resources = self._psij.ResourceSpecV1(process_count=mpi.process_count)
        return tuple(self.spec.mpi_argv), resources, self._mpi_launcher()

    def _project_metadata(self):
        if self._run_id is None:
            return None, None
        environment = {"ASE2SPRKKR_RUN_ID": self._run_id}
        attributes = None
        if self.executor.backend.lower() == "slurm":
            attributes = self._psij.JobAttributes(
                custom_attributes={"slurm.comment": f"ase2sprkkr:{self._run_id}"}
            )
        return environment, attributes

    def _create_job(self):
        argv, resources, launcher = self._job_command()
        if not argv:
            raise ValueError("Cannot submit a PSI/J run without a command")
        environment, attributes = self._project_metadata()
        job_spec = self._psij.JobSpec(
            executable=os.fspath(argv[0]),
            arguments=[os.fspath(arg) for arg in argv[1:]],
            directory=self.spec.directory,
            environment=environment,
            stdin_path=(
                self._stdin_path
                if not (self.spec.launch and self.spec.mpi_argv is not None)
                or self.spec.mpi_uses_stdin
                else None
            ),
            stdout_path=self._stdout_path,
            stderr_path=self._stderr_path,
            resources=resources,
            attributes=attributes,
            launcher=launcher,
        )
        return self._psij.Job(job_spec)

    @property
    def state(self):
        if self._exception is not None:
            return RunState.FAILED
        if self._observed_state is not None:
            return self._observed_state
        return _psij_state(self._job.status, attached=self._attached,
                           return_code=self.return_code)

    @property
    def exception(self):
        return self._exception

    @property
    def return_code(self):
        if self._observed_state is not None:
            # A completed wait() is authoritative even without an exit code.
            return self._return_code
        if self._return_code is not None:
            return self._return_code
        return getattr(getattr(self._job, "status", None), "exit_code", None)

    @property
    def failure_phase(self):
        return self._failure_phase

    @property
    def failure_message(self):
        return self._failure_message

    def set_manifest(self, files):
        self._files = files

    @property
    def files(self):
        return self._files

    def _update_files(self):
        if self._files is None:
            return
        self._files = discover_local_files(self._files)

    @property
    def backend_id(self):
        native_id = getattr(self._job, "native_id", None)
        return None if native_id is None else str(native_id)

    def add_output_listener(self, callback):
        self._output_observers.add(callback)

    def add_output_error_listener(self, callback):
        self._output_observers.add_error_listener(callback)

    def _prepare_io(self):
        if self.spec.input_file is not None and hasattr(self.spec.input_file, "flush"):
            self.spec.input_file.flush()
        protected = [
            _path(value, self.spec.directory)
            for value in (self.spec.input_file, self.spec.stdin, *self.spec.input_files,
                          *self.spec.cache_files)
        ]
        declared = [
            *(_path(value, self.spec.directory) for value in self.spec.output_files),
            self._stdout_path, self._stderr_path,
        ]
        validate_output_collisions(declared, protected=protected)
        for stream in (self.spec.output_file, self.spec.error_file):
            if stream is not None and hasattr(stream, "flush"):
                stream.flush()
        close_output_files(self.spec.output_file, self.spec.error_file)
        clear_declared_outputs(declared, protected=protected)
        for path in {self._stdout_path, self._stderr_path} - {None}:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_bytes(b"")

    async def _monitor(self):
        paths = {}
        for kind, path in (("out", self._stdout_path),
                           ("err", self._stderr_path)):
            # Merged stderr and stdout denote one physical stream.
            if path is not None and path not in paths.values():
                paths[kind] = path

        # PSI/J calls status callbacks on backend threads. Forward terminal
        # facts and explicit detach requests into this monitor's event loop.
        loop = asyncio.get_running_loop()
        terminal = loop.create_future()
        completed = []

        def publish(status):
            if not terminal.done():
                terminal.set_result(status)

        def schedule(status):
            try:
                if not loop.is_closed():
                    loop.call_soon_threadsafe(publish, status)
            except RuntimeError:
                # Loop shutdown after an observer has detached.
                pass

        def receive(status):
            if self._observation_stopped.is_set():
                return
            if _psij_state(status) not in (
                RunState.FINISHED, RunState.FAILED, RunState.CANCELLED,
            ):
                return
            completed.append(status)
            schedule(status)

        with self._observation_signal_lock:
            self._observation_wakeup = lambda: schedule(_OBSERVATION_STOPPED)
            if self._observation_stopped.is_set():
                self._observation_wakeup()

        detach = None
        read_error = None
        try:
            if not self._observation_stopped.is_set():
                detach = self.executor._observe_job(self._job, receive)
                # A job may finish before registration or during it.
                receive(self._job.status)
            if not paths:
                status = await terminal
            else:
                async def terminal_status():
                    return terminal.result() if terminal.done() else None

                status = await follow_growing_files(
                    paths, lambda kind, data: self._output_observers.emit_bytes(data, kind),
                    terminal_status, poll_interval=max(self.executor.poll_interval, 1.0),
                )
        except OSError as exc:
            read_error = exc
            status = completed[-1] if completed else None
        finally:
            if detach is not None:
                detach()
            with self._observation_signal_lock:
                self._observation_wakeup = None

        if status is _OBSERVATION_STOPPED:
            status = None

        observed_state = _psij_state(status) if status is not None else RunState.UNKNOWN
        if observed_state not in (
            RunState.FINISHED, RunState.FAILED, RunState.CANCELLED,
        ):
            self._observed_state = RunState.UNKNOWN
            self._failure_phase = None
            self._failure_message = (
                str(read_error) if read_error is not None
                else "PSI/J observation has no terminal job status"
            )
            # The job might still be writing its outputs. Do not finalize
            # their manifest on an inconclusive scheduler observation.
            return self.result

        self._return_code = getattr(status, "exit_code", None)
        # The scheduler may know FAILED/CANCELED without knowing the process
        # exit status; retain None rather than inventing an exit code.
        self.result.return_code = self._return_code
        self._observed_state = observed_state
        if observed_state == RunState.FAILED:
            self._failure_phase = FailurePhase.EXECUTION
            self._failure_message = getattr(status, "message", None) or (
                f"PSI/J job failed with exit code {self._return_code}"
                if self._return_code is not None else "PSI/J job failed"
            )
        if read_error is not None and observed_state == RunState.FINISHED:
            self._observed_state = RunState.FAILED
            self._failure_phase = FailurePhase.ARTIFACTS
            self._failure_message = str(read_error)
            return self.result
        if read_error is not None:
            # A failed or cancelled scheduler job remains terminal regardless
            # of whether its output stream could be read.
            return self.result
        try:
            self._update_files()
        except OSError as exc:
            if observed_state == RunState.FINISHED:
                self._observed_state = RunState.FAILED
                self._failure_phase = FailurePhase.ARTIFACTS
                self._failure_message = str(exc)
        return self.result

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
            self._failure_phase = self._failure_phase or FailurePhase.EXECUTION
            self._failure_message = self._failure_message or str(exc)
        finally:
            self._output_observers.finish_bytes()
            finish_monitor(self)

    def start(self):
        if self._attached:
            return start_monitor(self)
        return self.submit()

    def submit(self):
        with self._submit_lock:
            if self._attached:
                raise RuntimeError("Cannot submit an attached PSI/J job")
            if self._thread is not None or self.state != RunState.NEW:
                raise RuntimeError(f"Run is already {self.state.value}")
            try:
                self._prepare_io()
            except Exception as exc:
                # No scheduler interaction has happened yet.
                self._failure_phase = FailurePhase.SUBMISSION
                self._failure_message = str(exc)
                raise SubmissionNotStarted(str(exc)) from exc
            # An exception from PSI/J submit does not prove that the scheduler
            # rejected the job. Leave project submission recoverable even when
            # the native ID has not reached this process.
            self.executor._job_executor.submit(self._job)
            self.executor._cache_submitted_job(self._job)
            start_monitor(self)

    def wait(self):
        if self._thread is None:
            raise RuntimeError("Run has not been started")
        return wait_for_monitor(self)

    def run(self):
        if self._thread is None:
            self.start()
        return self.wait()

    async def run_async(self):
        if self._thread is None:
            self.start()
        return await asyncio.to_thread(self.wait)

    def stop_observation(self):
        """Stop monitoring without cancelling the scheduler calculation."""
        self._observation_stopped.set()
        with self._observation_signal_lock:
            wakeup = self._observation_wakeup
        if wakeup is not None:
            wakeup()

    def cancel(self):
        if self._thread is None and not self._attached:
            raise RuntimeError("Run has not been submitted")
        return request_cancellation(self, self._job.cancel)


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
        self._jobs = OrderedDict()
        self._observers = {}
        self._jobs_lock = RLock()
        self._max_idle_jobs = 64

    def create_run(self, spec: ExecutionSpec, callback=None):
        return PsijRunner(self, spec, callback=callback)

    def create_project_run(
        self, spec: ExecutionSpec, *, run_id, callback=None
    ):
        return PsijRunner(
            self,
            spec,
            callback=callback,
            run_id=run_id,
        )

    def status(self, *runs):
        from .endpoint import endpoint_for_executor

        expected = endpoint_for_executor(self)
        statuses = []
        for run in runs:
            if run.execution.endpoint != expected:
                raise ValueError("Run does not belong to this PSI/J endpoint")
            statuses.append(self.status_native(run.execution.backend_id))
        return tuple(statuses)

    def _trim_jobs(self):
        idle = sum(job not in self._observers for job in self._jobs.values())
        for native_id, job in tuple(self._jobs.items()):
            if idle <= self._max_idle_jobs:
                break
            if job not in self._observers:
                del self._jobs[native_id]
                idle -= 1

    def _cache_submitted_job(self, job):
        native_id = job.native_id
        if native_id is not None:
            with self._jobs_lock:
                self._jobs[str(native_id)] = job
                self._jobs.move_to_end(str(native_id))
                self._trim_jobs()

    def attach_native_job(self, native_id):
        native_id = str(native_id)
        with self._jobs_lock:
            job = self._jobs.get(native_id)
            if job is None:
                job = self._psij.Job()
                self._job_executor.attach(job, native_id)
                self._jobs[native_id] = job
            self._jobs.move_to_end(native_id)
            self._trim_jobs()
            return job

    def _observe_job(self, job, observer):
        with self._jobs_lock:
            observers = self._observers.get(job)
            if observers is None:
                # Register the first observer *before* installing the native
                # callback: PSI/J may call it synchronously. Roll back when
                # registration fails so a later attach can retry.
                observers = {observer}
                self._observers[job] = observers
                try:
                    job.set_job_status_callback(self._dispatch_job_status)
                except BaseException:
                    del self._observers[job]
                    raise
            else:
                observers.add(observer)

        def detach():
            with self._jobs_lock:
                current = self._observers.get(job)
                if current is None:
                    return
                current.discard(observer)
                if current:
                    return
                del self._observers[job]
                # PSI/J 0.9.11 incorrectly wraps None as a callback.
                # This is best-effort teardown: failing to unregister must
                # not replace a monitor's original failure or its terminal
                # result. The dispatcher has no subscribers after removal.
                try:
                    job.set_job_status_callback(_ignore_job_status)
                except Exception:
                    _log.exception("Failed to disable PSI/J job status callback")
                native_id = job.native_id
                if (native_id is not None
                        and self._jobs.get(str(native_id)) is job
                        and _psij_state(job.status) in (
                            RunState.FINISHED, RunState.FAILED, RunState.CANCELLED,
                        )):
                    del self._jobs[str(native_id)]
                self._trim_jobs()
        return detach

    def _dispatch_job_status(self, job, status):
        with self._jobs_lock:
            observers = tuple(self._observers.get(job, ()))
        for observer in observers:
            observer(status)

    def _status_job(self, native_id):
        return self.attach_native_job(native_id)

    def status_native(self, native_id):
        """Return the currently known native job status without waiting."""
        native_id = str(native_id)
        if self.backend.lower() == "slurm" and self.url is None:
            return _slurm_status_native(native_id)

        job = self._status_job(native_id)
        status = getattr(job, "status", None)
        state = _psij_state(status, attached=True)
        return_code = getattr(status, "exit_code", None)
        if state in (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED):
            with self._jobs_lock:
                if job not in self._observers and self._jobs.get(native_id) is job:
                    del self._jobs[native_id]
        return RunStatus(
            state=state, backend_id=native_id,
            return_code=return_code,
            message=getattr(status, "message", None),
            failure_phase=FailurePhase.EXECUTION if state == RunState.FAILED else None,
        )

    def wait_status_native(self, native_id, timeout):
        """Wait up to *timeout* seconds for a known state of a native job.

        Unlike status_native(), this method may wait. A timeout returns UNKNOWN
        if the scheduler has still not provided a meaningful job state.
        """
        timeout = float(timeout)
        if timeout < 0:
            raise ValueError("Status timeout must not be negative")
        status = self.status_native(native_id)
        if status.state != RunState.UNKNOWN or timeout == 0:
            return status

        job = self._status_job(native_id)
        # The PSI/J Job object waits for notifications from its executor.
        # There is no repeated scheduler status query or sleep loop.
        job.wait(
            timeout=timedelta(seconds=timeout),
            target_states=[self._psij.JobState.QUEUED],
        )
        return self.status_native(native_id)

    def recover_native(self, run_id, directory, *, submitted_at=None):
        """Return Slurm's native recovery facts without a client RunRecord."""
        if self.backend.lower() != "slurm" or self.url is not None:
            raise NotImplementedError(
                f"PSI/J backend {self.backend!r} at {self.url!r} cannot "
                "recover uncertain submissions using local Slurm commands"
            )
        matches = _slurm_current_matches(run_id)
        if len(matches) > 1:
            raise RuntimeError(
                f"Multiple live Slurm jobs carry ASE2SPRKKR run id {run_id!r}"
            )
        if matches:
            native_id, state, workdir = matches[0]
            try:
                details = _slurm_details(native_id)
            except (FileNotFoundError, subprocess.CalledProcessError):
                # The job may have finished between squeue and scontrol.
                # Recover its terminal facts from accounting instead.
                details = None
            if details:
                return {
                    "native_id": str(native_id), "state": _slurm_state(state),
                    "directory": details.get("WorkDir") or workdir or directory,
                    "stdout": details.get("StdOut"),
                    "stderr": details.get("StdErr"), "return_code": None,
                }

        if submitted_at is None or not _slurm_accounting_stores_comments():
            raise NotImplementedError(
                "Slurm cannot safely determine whether this submission already "
                "finished because job comments are unavailable in accounting"
            )
        accounting = _slurm_accounting_matches(run_id, submitted_at)
        if len(accounting) > 1:
            raise RuntimeError(
                f"Multiple Slurm accounting records carry ASE2SPRKKR run id {run_id!r}"
            )
        if not accounting:
            # An empty sacct result is not authoritative negative evidence:
            # Slurm accounting can lag behind the scheduler. Do not release
            # the client's SUBMITTING reservation and risk a duplicate job.
            raise NotImplementedError(
                "Slurm has no matching job in squeue or accounting yet; "
                "submission acceptance is still uncertain"
            )
        native_id, state, workdir, stdout, stderr, return_code = accounting[0]
        recovered_state = _slurm_state(state)
        if recovered_state == RunState.FINISHED and return_code not in (None, 0):
            recovered_state = RunState.FAILED
        return {
            "native_id": str(native_id), "state": recovered_state,
            "directory": workdir or directory, "stdout": stdout,
            "stderr": stderr, "return_code": return_code,
        }

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        facts = self.recover_native(
            run_id, submission.execution.directory.target, submitted_at=submitted_at,
        )
        if facts is None:
            return None
        execution = submission.execution.submitted(facts["native_id"])
        directory = facts["directory"]
        if directory != execution.directory.target:
            from .run_model import WorkingDirectory
            execution = replace(
                execution,
                directory=WorkingDirectory(execution.directory.local, directory),
            )
        files = submission.files.with_remote_stream_paths(
            stdout=facts["stdout"], stderr=facts["stderr"],
        )
        state = facts["state"]
        return RunStatus(
            state=state, backend_id=facts["native_id"],
            execution=execution, files=files,
            return_code=facts["return_code"],
            failure_phase=FailurePhase.EXECUTION if state == RunState.FAILED else None,
        )

    def attach(self, run, callback=None, output_listener=None, start=True):
        from .endpoint import endpoint_for_executor

        if run.execution.endpoint != endpoint_for_executor(self):
            raise ValueError("Run does not belong to this PSI/J endpoint")
        backend_id = run.execution.backend_id
        if not backend_id:
            raise ValueError("PSI/J run lacks backend ID")
        job = self.attach_native_job(backend_id)
        spec = ExecutionSpec(
            argv=(),
            input_file=run.files.input_file,
            output_file=run.files.output_file,
            error_file=None if run.files.stderr is None else run.files.stderr.local_path,
            directory=run.execution.directory.local,
        )
        runner = PsijRunner(
            self, spec, callback=callback, output_listener=output_listener,
            start=False, job=job, attached=True,
        )
        runner.set_manifest(run.files)
        if start:
            runner.start()
        return runner
