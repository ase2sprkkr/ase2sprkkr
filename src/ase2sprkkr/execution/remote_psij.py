"""PSI/J backend used by the detached remote execution helper."""

from __future__ import annotations

import asyncio
import json
import os
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from .base import ExecutionSpec, LaunchSpec, RunReference, RunState
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
    """Return an exclusively locked file, or ``None`` if an agent owns it."""
    import fcntl

    path = _agent_lock_path(directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.open("a+b")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        return None
    return lock


def _agent_is_running(directory):
    import fcntl

    lock = _try_agent_lock(directory)
    if lock is None:
        return True
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    finally:
        lock.close()
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
    """Ensure an agent exists when any requested PSI/J run is still active."""
    from .remote import _read_state, _request_path

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
    """Own shared PSI/J executors and all active PSI/J runs in one workspace."""

    def __init__(self, directory):
        self.directory = Path(directory).expanduser().resolve()
        self._executors = {}
        self._active = {}

    def _executor_for(self, request):
        key = _executor_key(request)
        executor = self._executors.get(key)
        if executor is None:
            executor = _executor(request)
            self._executors[key] = executor
        return executor

    def _runs(self):
        from .remote import _request_path, list_remote_runs

        for state in list_remote_runs(self.directory):
            if state.get("state") in _FINAL_STATE_NAMES:
                continue
            run_id = state.get("id")
            if not run_id:
                continue
            try:
                request = json.loads(_request_path(self.directory, run_id).read_text())
            except (FileNotFoundError, json.JSONDecodeError):
                continue
            if request.get("executor", "local") == "psij":
                yield run_id, state, request

    def _fail(self, run_id, exc):
        from .remote import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        with _state_lock(self.directory, run_id):
            state = _read_state_unlocked(self.directory, run_id)
            if state.get("state") not in _FINAL_STATE_NAMES:
                _write_state_unlocked(
                    self.directory,
                    run_id,
                    state=RunState.FAILED.value,
                    error=repr(exc),
                    return_code=None,
                    finished_at=state.get("finished_at") or _timestamp(),
                )
        self._active.pop(run_id, None)

    def _submit(self, run_id, request):
        from .remote import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        executor = self._executor_for(request)
        output_file = request.get("output_file") or os.path.join(self.directory, "stdout")
        error_file = request.get("error_file") or output_file + ".err"
        spec = _execution_spec(self.directory, request, output_file, error_file)
        runner = executor.create_run(_RemoteReader(), spec)

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
            if state.get("state") != RunState.QUEUED.value:
                raise RuntimeError(
                    f"Cannot submit remote PSI/J run {run_id!r} from state "
                    f"{state.get('state')!r}"
                )

            # Keep the per-run state lock across scheduler submission. This is
            # the handoff boundary with cancellation: a cancel request either
            # wins before submit, or observes a fully published native id.
            runner.submit()
            reference = runner.reference
            if reference is None:
                raise RuntimeError("PSI/J did not provide a persistent job reference")

            current = runner.state
            published = (
                current
                if current in (RunState.QUEUED, RunState.RUNNING)
                else RunState.QUEUED
            )
            _write_state_unlocked(
                self.directory,
                run_id,
                state=published.value,
                pid=os.getpid(),
                backend_reference=asdict(reference),
                started_at=state.get("started_at") or _timestamp(),
            )

        self._active[run_id] = _ActiveRun(runner, current)
        return True

    def _attach(self, run_id, state, request):
        reference = state.get("backend_reference")
        if reference is None:
            return False
        executor = self._executor_for(request)
        runner = executor.attach(RunReference(**reference), _RemoteReader())
        self._active[run_id] = _ActiveRun(
            runner,
            RunState(state.get("state", RunState.UNKNOWN.value)),
        )
        return True

    def _discover(self):
        found = False
        for run_id, state, request in self._runs():
            found = True
            if run_id in self._active:
                continue
            try:
                if state.get("backend_reference") is not None:
                    self._attach(run_id, state, request)
                elif state.get("state") == RunState.QUEUED.value:
                    self._submit(run_id, request)
            except BaseException as exc:
                self._fail(run_id, exc)
        return found

    def _finish(self, run_id, active, requested_cancel):
        from .remote import (
            _read_state_unlocked,
            _state_lock,
            _timestamp,
            _write_state_unlocked,
        )

        try:
            result = active.runner.wait()
            return_code = result.return_code
            final = active.runner.state
            if final not in _FINAL_RUN_STATES:
                final = RunState.FINISHED if return_code == 0 else RunState.FAILED
        except BaseException as exc:
            if requested_cancel:
                final = RunState.CANCELLED
                return_code = None
            else:
                self._fail(run_id, exc)
                return

        with _state_lock(self.directory, run_id):
            state = _read_state_unlocked(self.directory, run_id)
            if state.get("cancel_requested") or requested_cancel:
                final = RunState.CANCELLED
            _write_state_unlocked(
                self.directory,
                run_id,
                state=final.value,
                return_code=return_code,
                finished_at=state.get("finished_at") or _timestamp(),
            )
        self._active.pop(run_id, None)

    def _poll_active(self):
        from .remote import (
            _read_state,
            _read_state_unlocked,
            _state_lock,
            _write_state_unlocked,
        )

        for run_id, active in tuple(self._active.items()):
            try:
                state = _read_state(self.directory, run_id)
                requested_cancel = bool(state.get("cancel_requested"))
                current = active.runner.state
                if current in _FINAL_RUN_STATES:
                    self._finish(run_id, active, requested_cancel)
                    continue

                if requested_cancel and not active.cancel_sent:
                    active.runner.cancel()
                    active.cancel_sent = True
                    current = active.runner.state
                    if current in _FINAL_RUN_STATES:
                        self._finish(run_id, active, requested_cancel)
                        continue

                if current in (RunState.QUEUED, RunState.RUNNING) and current != active.last_state:
                    with _state_lock(self.directory, run_id):
                        current_state = _read_state_unlocked(self.directory, run_id)
                        if current_state.get("state") not in _FINAL_STATE_NAMES:
                            _write_state_unlocked(
                                self.directory, run_id, state=current.value
                            )
                active.last_state = current
            except BaseException as exc:
                self._fail(run_id, exc)

    def run(self):
        idle_since = None
        while True:
            found = self._discover()
            self._poll_active()
            if found or self._active:
                idle_since = None
            elif idle_since is None:
                idle_since = time.monotonic()
            elif time.monotonic() - idle_since >= _AGENT_IDLE_TIMEOUT:
                return
            time.sleep(_AGENT_POLL_INTERVAL)


def run_remote_psij_worker(directory, run_id, request):
    """Join or become the workspace PSI/J agent.

    The first detached worker owns the workspace agent lock and stays alive to
    manage all PSI/J runs. Later per-submit workers see the lock, exit quickly,
    and leave scheduler polling to the shared PSI/J executor instance.
    """
    directory = Path(directory).expanduser().resolve()
    lock = _try_agent_lock(directory)
    if lock is None:
        return
    try:
        _RemotePsijAgent(directory).run()
    finally:
        import fcntl

        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()


def cancel_remote_psij(directory, run_id, request, timeout=1.0):
    """Request cancellation through the shared target-side PSI/J agent."""
    from .remote import (
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
