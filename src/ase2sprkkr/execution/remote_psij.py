"""PSI/J backend used by the detached execution helper."""

from __future__ import annotations

import asyncio
import json
import os
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace

from ..common.file_locks import FileLock
from .base import ExecutionSpec, FailurePhase, LaunchSpec, RunReference, RunState
from .psij import PsijExecutor


_FINAL_RUN_STATES = {
    RunState.FINISHED,
    RunState.FAILED,
    RunState.CANCELLED,
}
_FINAL_STATE_NAMES = {state.value for state in _FINAL_RUN_STATES}
_AGENT_LOCK_FILE = "psij-agent.lock"
_AGENT_POLL_INTERVAL = 0.1
_AGENT_IDLE_TIMEOUT = 5.0
_RETRY_INITIAL = 0.5
_RETRY_MAX = 10.0


class _RemoteResult:
    return_code = None


class _RemoteReader:
    """Drain scheduler output; parsing is done by the local RemoteRunner."""

    print_output = False

    def create_result(self, output_file, input_file):
        return _RemoteResult()

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.return_code = return_code
        return result


def _executor(request):
    options = request.get("executor_options") or {}
    return PsijExecutor(
        options.get("backend") or "slurm",
        url=options.get("url"),
    )


def _executor_key(request):
    """Stable key for one target-side PSI/J executor configuration."""
    return json.dumps(
        request.get("executor_options") or {},
        sort_keys=True,
        separators=(",", ":"),
        default=repr,
    )


def _execution_spec(directory, request, output_file, error_file):
    mpi_argv = request.get("mpi_argv")
    return ExecutionSpec(
        argv=tuple(request["argv"]),
        input_file=request.get("input_file"),
        output_file=output_file,
        error_file=error_file,
        directory=os.fspath(directory),
        stdin=request.get("stdin"),
        mpi_argv=None if mpi_argv is None else tuple(mpi_argv),
        launch=LaunchSpec(**(request.get("launch") or {})),
    )


def _agent_lock_path(directory):
    return Path(directory) / ".ase2sprkkr" / _AGENT_LOCK_FILE


def _try_agent_lock(directory):
    """Return an acquired lock, or ``None`` if an agent owns it."""
    lock = FileLock(_agent_lock_path(directory))
    if not lock.acquire(blocking=False):
        return None
    return lock


def _agent_is_running(directory):
    lock = _try_agent_lock(directory)
    if lock is None:
        return True
    lock.release()
    return False


def _spawn_agent(directory, run_id):
    subprocess.Popen(
        [
            sys.executable,
            "-m",
            "ase2sprkkr.execution.remote_worker",
            os.fspath(directory),
            str(run_id),
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
        close_fds=True,
    )


def ensure_remote_psij_agent(directory, run_ids):
    """Ensure an agent exists when a requested PSI/J run needs attention."""
    from .detached import _read_state, _request_path

    directory = Path(directory).expanduser().resolve()
    if _agent_is_running(directory):
        return True

    for run_id in run_ids:
        try:
            state = _read_state(directory, run_id)
            request = json.loads(_request_path(directory, run_id).read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            continue
        if state.get("state") in _FINAL_STATE_NAMES:
            continue
        if request.get("executor", "local") != "psij":
            continue
        _spawn_agent(directory, run_id)
        return True
    return False


@dataclass
class _ActiveRun:
    runner: object
    last_state: RunState
    cancel_sent: bool = False


class _RemotePsijAgent:
    """Temporarily own shared PSI/J executors for one detached workspace."""

    def __init__(self, directory):
        self.directory = Path(directory).expanduser().resolve()
        self._executors = {}
        self._active = {}
        self._retry = {}

    def _retry_ready(self, run_id):
        retry = self._retry.get(run_id)
        return retry is None or time.monotonic() >= retry[1]

    def _schedule_retry(self, run_id):
        attempts = self._retry.get(run_id, (0, 0.0))[0] + 1
        exponent = min(attempts - 1, 10)
        delay = min(_RETRY_INITIAL * (2**exponent), _RETRY_MAX)
        self._retry[run_id] = (attempts, time.monotonic() + delay)

    def _clear_retry(self, run_id):
        self._retry.pop(run_id, None)

    def _executor_for(self, request):
        key = _executor_key(request)
        executor = self._executors.get(key)
        if executor is None:
            executor = _executor(request)
            self._executors[key] = executor
        return executor

    def _runs(self):
        from .detached import _request_path, list_detached_runs

        current = set()
        for state in list_detached_runs(self.directory):
            run_id = state.get("id")
            if run_id:
                current.add(run_id)
            if state.get("state") in _FINAL_STATE_NAMES:
                if run_id:
                    self._clear_retry(run_id)
                continue
            if not run_id:
                continue
            try:
                request = json.loads(_request_path(self.directory, run_id).read_text())
            except (FileNotFoundError, json.JSONDecodeError):
                continue
            if request.get("executor", "local") == "psij":
                yield run_id, state, request
        for run_id in tuple(self._retry):
            if run_id not in current:
                self._clear_retry(run_id)

    def _fail(self, run_id, exc):
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        with _state_lock(self.directory, run_id):
            state = _read_state_unlocked(self.directory, run_id)
            if state.get("state") not in _FINAL_STATE_NAMES:
                phase = (
                    FailurePhase.EXECUTION
                    if state.get("started_at") is not None
                    else FailurePhase.SUBMISSION
                    if state.get("backend_reference") is None
                    else FailurePhase.LAUNCH
                )
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    state=RunState.FAILED.value,
                    error=repr(exc),
                    failure_phase=phase.value,
                    return_code=None,
                    finished_at=state.get("finished_at") or _timestamp(),
                )
        self._active.pop(run_id, None)
        self._clear_retry(run_id)

    def _record_retryable_error(self, run_id, field, exc):
        """Persist retry diagnostics when possible without changing lifecycle."""
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _write_state_unlocked,
        )

        try:
            with _state_lock(self.directory, run_id):
                state = _read_state_unlocked(self.directory, run_id)
                if state.get("state") not in _FINAL_STATE_NAMES:
                    _write_state_unlocked(
                        self.directory,
                        run_id,
                        **{field: repr(exc), "error": repr(exc)},
                    )
        except OSError:
            return False
        return True

    def _clear_retryable_error(self, run_id, field):
        """Clear one recovered diagnostic without overwriting unrelated errors."""
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _write_state_unlocked,
        )

        try:
            with _state_lock(self.directory, run_id):
                state = _read_state_unlocked(self.directory, run_id)
                if state.get("state") in _FINAL_STATE_NAMES:
                    return True
                message = state.get(field)
                if message is None:
                    return True
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    **{
                        field: None,
                        "error": None if state.get("error") == message else state.get("error"),
                    },
                )
        except OSError:
            return False
        return True

    def _submit(self, run_id, request):
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        executor = self._executor_for(request)
        output_file = request.get("output_file") or os.path.join(self.directory, "stdout")
        error_file = request.get("error_file") or output_file + ".err"
        spec = _execution_spec(self.directory, request, output_file, error_file)
        runner = executor.create_persistent_run(_RemoteReader(), spec, run_id=run_id)

        with _state_lock(self.directory, run_id):
            state = _read_state_unlocked(self.directory, run_id)
            if state.get("state") in _FINAL_STATE_NAMES:
                return False
            if state.get("backend_reference") is not None:
                return False
            if state.get("cancel_requested"):
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    state=RunState.CANCELLED.value,
                    finished_at=state.get("finished_at") or _timestamp(),
                )
                return True
            if (
                state.get("state") != RunState.SUBMITTING.value
                or state.get("submitted_at") is not None
            ):
                raise RuntimeError(
                    f"Cannot submit remote PSI/J run {run_id!r} from state "
                    f"{state.get('state')!r} with submitted_at="
                    f"{state.get('submitted_at')!r}"
                )
            submitted_at = _timestamp()
            _write_state_unlocked(
                self.directory,
                run_id,
                pid=os.getpid(),
                submitted_at=submitted_at,
                started_at=None,
                submission_error=None,
                recovery_error=None,
                error=None,
                failure_phase=None,
            )

        try:
            runner.submit()
        except BaseException as exc:
            # Submission failure can be ambiguous: the scheduler may already
            # have accepted the job even if the client call raised. Keep the
            # explicit SUBMITTING state and let recovery by persistent run id
            # decide whether a scheduler job exists before releasing the run.
            self._record_retryable_error(run_id, "submission_error", exc)
            self._schedule_retry(run_id)
            return True

        try:
            reference = runner.reference
            if reference is None:
                exc = RuntimeError(
                    "PSI/J did not provide a persistent job reference"
                )
                self._record_retryable_error(run_id, "recovery_error", exc)
                self._schedule_retry(run_id)
                return True

            current = runner.state
            published = (
                current
                if current in (RunState.QUEUED, RunState.RUNNING)
                else RunState.QUEUED
            )
            with _state_lock(self.directory, run_id):
                latest = _read_state_unlocked(self.directory, run_id)
                now = _timestamp()
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    state=published.value,
                    pid=os.getpid(),
                    backend_reference=asdict(reference),
                    submitted_at=latest.get("submitted_at") or now,
                    started_at=(
                        latest.get("started_at") or now
                        if published == RunState.RUNNING
                        else latest.get("started_at")
                    ),
                    submission_error=None,
                    recovery_error=None,
                    attach_error=None,
                    error=None,
                    failure_phase=None,
                )
        except BaseException as exc:
            # Once submit() returned, a scheduler job may already exist. Any
            # failure while reading or persisting its identity is therefore an
            # uncertain handoff, not a terminal submission/launch failure.
            # Leave SUBMITTING intact and recover by persistent run id.
            self._record_retryable_error(run_id, "recovery_error", exc)
            self._schedule_retry(run_id)
            return True

        self._clear_retry(run_id)
        self._active[run_id] = _ActiveRun(runner, current)
        return True

    def _recover_submitting(self, run_id, state, request):
        """Recover the scheduler identity after a target-agent submit crash."""
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        submission = SimpleNamespace(
            directory=request.get("directory") or os.fspath(self.directory)
        )
        try:
            executor = self._executor_for(request)
            status = executor.recover_submission(
                run_id,
                submission,
                submitted_at=state.get("submitted_at"),
            )
        except NotImplementedError:
            # Absence cannot be proven. Leaving SUBMITTING is deliberately safer
            # than creating a duplicate scheduler job.
            return False
        except Exception as exc:
            self._record_retryable_error(run_id, "recovery_error", exc)
            self._schedule_retry(run_id)
            return False

        if status is None:
            self._clear_retry(run_id)
            with _state_lock(self.directory, run_id):
                latest = _read_state_unlocked(self.directory, run_id)
                if (
                    latest.get("state") != RunState.SUBMITTING.value
                    or latest.get("backend_reference") is not None
                ):
                    return False
                if latest.get("cancel_requested"):
                    _write_state_unlocked(
                        self.directory,
                        run_id,
                        state=RunState.CANCELLED.value,
                        finished_at=latest.get("finished_at") or _timestamp(),
                        recovery_error=None,
                    )
                elif latest.get("submission_error"):
                    _write_state_unlocked(
                        self.directory,
                        run_id,
                        state=RunState.FAILED.value,
                        failure_phase=FailurePhase.SUBMISSION.value,
                        return_code=None,
                        finished_at=latest.get("finished_at") or _timestamp(),
                        recovery_error=None,
                    )
                else:
                    # A post-submit handoff failed, but recovery proves that no
                    # scheduler job exists. Reset only the submission attempt;
                    # the next discovery pass may submit it again.
                    _write_state_unlocked(
                        self.directory,
                        run_id,
                        state=RunState.SUBMITTING.value,
                        pid=None,
                        submitted_at=None,
                        recovery_error=None,
                        error=None,
                    )
            return True

        try:
            reference = status.reference
            if reference is None:
                raise RuntimeError("Recovered PSI/J submission lacks a job reference")
            recovered_state = RunState(status.state)
            phase = getattr(status, "failure_phase", None)
            execution_evidence = recovered_state in (
                RunState.RUNNING,
                RunState.FINISHED,
            )
            with _state_lock(self.directory, run_id):
                latest = _read_state_unlocked(self.directory, run_id)
                if latest.get("state") != RunState.SUBMITTING.value:
                    return False
                now = _timestamp()
                changes = {
                    "state": recovered_state.value,
                    "backend_reference": asdict(reference),
                    "pid": os.getpid(),
                    "submission_error": None,
                    "recovery_error": None,
                    "error": status.message,
                    "failure_phase": None if phase is None else phase.value,
                    "submitted_at": latest.get("submitted_at") or now,
                }
                if execution_evidence:
                    changes["started_at"] = latest.get("started_at") or now
                if status.return_code is not None:
                    changes["return_code"] = status.return_code
                if recovered_state in _FINAL_RUN_STATES:
                    if recovered_state == RunState.FAILED and phase is None:
                        changes["failure_phase"] = (
                            FailurePhase.EXECUTION.value
                            if execution_evidence
                            else FailurePhase.LAUNCH.value
                        )
                    changes["finished_at"] = latest.get("finished_at") or now
                recovered = _write_state_unlocked(self.directory, run_id, **changes)
        except Exception as exc:
            # A non-None recovery status proves that the scheduler knows this
            # job. Failure to normalize or persist the recovered identity must
            # therefore stay retryable; treating it as SUBMISSION failure could
            # release the persistent namespace and create a duplicate job.
            self._record_retryable_error(run_id, "recovery_error", exc)
            self._schedule_retry(run_id)
            return False

        self._clear_retry(run_id)
        if recovered_state in _FINAL_RUN_STATES:
            return True
        return self._attach(run_id, recovered, request)

    def _attach(self, run_id, state, request):
        reference = state.get("backend_reference")
        if reference is None:
            return False
        try:
            executor = self._executor_for(request)
            runner = executor.attach(RunReference(**reference), _RemoteReader())
        except Exception as exc:
            self._record_retryable_error(run_id, "attach_error", exc)
            self._schedule_retry(run_id)
            return True
        self._clear_retry(run_id)
        self._active[run_id] = _ActiveRun(
            runner,
            RunState(state.get("state", RunState.UNKNOWN.value)),
        )
        from .detached import _read_state_unlocked, _state_lock, _write_state_unlocked

        with _state_lock(self.directory, run_id):
            latest = _read_state_unlocked(self.directory, run_id)
            if latest.get("state") not in _FINAL_STATE_NAMES:
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    attach_error=None,
                    error=(
                        None
                        if latest.get("error") == latest.get("attach_error")
                        else latest.get("error")
                    ),
                )
        return True

    def _discover(self):
        activity = False
        for run_id, state, request in self._runs():
            if run_id in self._active or not self._retry_ready(run_id):
                continue
            try:
                if state.get("backend_reference") is not None:
                    activity = self._attach(run_id, state, request) or activity
                elif state.get("state") == RunState.SUBMITTING.value:
                    if state.get("submitted_at") is None:
                        activity = self._submit(run_id, request) or activity
                    else:
                        activity = (
                            self._recover_submitting(run_id, state, request) or activity
                        )
                elif state.get("state") == RunState.QUEUED.value:
                    raise RuntimeError(
                        f"Queued PSI/J run {run_id!r} has no backend reference"
                    )
            except BaseException as exc:
                self._fail(run_id, exc)
                activity = True
        return activity

    def _finish(self, run_id, active, requested_cancel):
        from .detached import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        try:
            result = active.runner.wait()
            return_code = getattr(active.runner, "return_code", None)
            if return_code is None:
                return_code = result.return_code
            final = active.runner.state
        except BaseException as exc:
            final = active.runner.state
            return_code = getattr(active.runner, "return_code", None)
            if requested_cancel:
                final = RunState.CANCELLED
            elif final not in _FINAL_RUN_STATES:
                self._fail(run_id, exc)
                return

        if final not in _FINAL_RUN_STATES:
            final = RunState.FINISHED if return_code == 0 else RunState.FAILED

        with _state_lock(self.directory, run_id):
            state = _read_state_unlocked(self.directory, run_id)
            if state.get("cancel_requested") or requested_cancel:
                final = RunState.CANCELLED
            now = _timestamp()
            execution_evidence = (
                state.get("started_at") is not None
                or active.last_state == RunState.RUNNING
                or final == RunState.FINISHED
            )
            changes = {
                "state": final.value,
                "return_code": return_code,
                "finished_at": state.get("finished_at") or now,
            }
            if execution_evidence:
                changes["started_at"] = state.get("started_at") or now
            if final == RunState.FAILED:
                changes["failure_phase"] = (
                    FailurePhase.EXECUTION.value
                    if execution_evidence
                    else FailurePhase.LAUNCH.value
                )
            elif final != RunState.CANCELLED:
                changes["failure_phase"] = None
            _write_state_unlocked(self.directory, run_id, **changes)
        self._active.pop(run_id, None)
        self._clear_retry(run_id)

    def _poll_active(self):
        from .detached import (
            _read_state,
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        activity = False
        for run_id, active in tuple(self._active.items()):
            if not self._retry_ready(run_id):
                continue

            try:
                state = _read_state(self.directory, run_id)
            except OSError as exc:
                self._record_retryable_error(run_id, "poll_error", exc)
                self._schedule_retry(run_id)
                activity = True
                continue

            requested_cancel = bool(state.get("cancel_requested"))
            try:
                current = active.runner.state
            except Exception as exc:
                self._record_retryable_error(run_id, "poll_error", exc)
                self._schedule_retry(run_id)
                activity = True
                continue

            if current in _FINAL_RUN_STATES:
                self._clear_retry(run_id)
                self._clear_retryable_error(run_id, "poll_error")
                self._finish(run_id, active, requested_cancel)
                activity = True
                continue

            if requested_cancel and not active.cancel_sent:
                try:
                    active.runner.cancel()
                    current = active.runner.state
                except Exception as exc:
                    self._record_retryable_error(run_id, "poll_error", exc)
                    self._schedule_retry(run_id)
                    activity = True
                    continue
                active.cancel_sent = True
                activity = True
                if current in _FINAL_RUN_STATES:
                    self._clear_retry(run_id)
                    self._clear_retryable_error(run_id, "poll_error")
                    self._finish(run_id, active, requested_cancel)
                    continue

            if (
                current in (RunState.QUEUED, RunState.RUNNING)
                and current != active.last_state
            ):
                try:
                    with _state_lock(self.directory, run_id):
                        current_state = _read_state_unlocked(self.directory, run_id)
                        if current_state.get("state") not in _FINAL_STATE_NAMES:
                            changes = {"state": current.value}
                            if current == RunState.RUNNING:
                                changes["started_at"] = (
                                    current_state.get("started_at") or _timestamp()
                                )
                            _write_state_unlocked(
                                self.directory, run_id, **changes
                            )
                except OSError as exc:
                    self._record_retryable_error(run_id, "poll_error", exc)
                    self._schedule_retry(run_id)
                    activity = True
                    continue
                activity = True

            active.last_state = current
            self._clear_retry(run_id)
            self._clear_retryable_error(run_id, "poll_error")
        return activity

    def run(self):
        idle_since = None
        while True:
            activity = self._discover()
            activity = self._poll_active() or activity
            if self._active or activity or self._retry:
                idle_since = None
            elif idle_since is None:
                idle_since = time.monotonic()
            elif time.monotonic() - idle_since >= _AGENT_IDLE_TIMEOUT:
                return
            time.sleep(_AGENT_POLL_INTERVAL)


def run_remote_psij_worker(directory, run_id, request):
    """Join or become the workspace PSI/J agent."""
    directory = Path(directory).expanduser().resolve()
    lock = _try_agent_lock(directory)
    if lock is None:
        return
    try:
        _RemotePsijAgent(directory).run()
    finally:
        lock.release()


def cancel_remote_psij(directory, run_id, request, timeout=1.0):
    """Request cancellation through the target-side PSI/J agent."""
    from .detached import (
        _read_state,
        _read_state_unlocked,
        _state_lock,
        _write_state_unlocked,
    )

    directory = Path(directory).expanduser().resolve()
    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATE_NAMES:
            return state
        _write_state_unlocked(
            directory,
            run_id,
            cancel_requested=True,
        )

    ensure_remote_psij_agent(directory, (run_id,))

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        state = _read_state(directory, run_id)
        if state.get("state") in _FINAL_STATE_NAMES:
            return state
        time.sleep(0.05)
    return _read_state(directory, run_id)
