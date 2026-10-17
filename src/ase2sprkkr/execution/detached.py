"""Crash-durable detached protocol for persistent local execution.

This protocol is intentionally local-only. Remote execution uses the stateless
``remote_service`` and never stores ASE2SPRKKR lifecycle records on the target.
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ..common.file_locks import FileLock
from ._durability import durable_mkdir, durable_write_json, sync_file
from .base import ExecutionSpec, FailurePhase, LaunchSpec, RunState


_FINAL_STATES = {"finished", "failed", "cancelled"}
_RESULT_FINAL_STATES = {"finished", "failed"}
_METADATA_DIR = ".ase2sprkkr"
_RUNS_DIR = "runs"
_STATE_FILE = "state.json"
_STATE_LOCK_FILE = "state.lock"
_REQUEST_FILE = "request.json"


def detached_run_state(value):
    if value == "failed":
        return RunState.FAILED
    return RunState(value)


def detached_failure_phase(state):
    value = state.get("failure_phase")
    return None if value is None else FailurePhase(value)


def detached_launch_failed(state):
    return (
        detached_run_state(state.get("state", "unknown")) == RunState.FAILED
        and detached_failure_phase(state) == FailurePhase.LAUNCH
    )


def detached_stream_file(state, key):
    stream = state.get(f"stream_{key}")
    if stream and os.path.exists(stream):
        return stream
    return state.get(key)


def _timestamp():
    return datetime.now(timezone.utc).isoformat()


def _run_directory(directory, run_id):
    return Path(directory) / _METADATA_DIR / _RUNS_DIR / run_id


def _state_path(directory, run_id):
    return _run_directory(directory, run_id) / _STATE_FILE


def _state_lock_path(directory, run_id):
    return _run_directory(directory, run_id) / _STATE_LOCK_FILE


def _request_path(directory, run_id):
    return _run_directory(directory, run_id) / _REQUEST_FILE


def _state_lock(directory, run_id):
    return FileLock(_state_lock_path(directory, run_id))


def _read_state_unlocked(directory, run_id):
    with _state_path(directory, run_id).open() as file:
        return json.load(file)


def _read_state(directory, run_id):
    return _read_state_unlocked(directory, run_id)


def _sync_result_files_before_final_state(state):
    """Make local terminal publication a durable filesystem boundary."""
    if state.get("state") not in _RESULT_FINAL_STATES:
        return
    if state.get("failure_phase") in {
        FailurePhase.SUBMISSION.value,
        FailurePhase.LAUNCH.value,
        FailurePhase.ARTIFACTS.value,
    }:
        return

    seen = set()
    for key in ("output_file", "error_file"):
        path = state.get(key)
        if path is None:
            continue
        path = os.fspath(path)
        if path in seen:
            continue
        seen.add(path)
        if os.path.isfile(path):
            sync_file(path)


def _write_state_unlocked(workspace, run_id, **changes):
    path = _state_path(workspace, run_id)
    try:
        state = _read_state_unlocked(workspace, run_id)
    except FileNotFoundError:
        state = {}
    state.update(changes)
    try:
        _sync_result_files_before_final_state(state)
    except OSError as exc:
        failure = dict(state)
        failure.update(
            state="failed",
            error=repr(exc),
            failure_phase=FailurePhase.ARTIFACTS.value,
            finished_at=state.get("finished_at") or _timestamp(),
        )
        durable_write_json(path, failure)
        raise
    durable_write_json(path, state)
    return state


def _write_state(workspace, run_id, **changes):
    with _state_lock(workspace, run_id):
        return _write_state_unlocked(workspace, run_id, **changes)


def _absolute_file(directory, path):
    if path is None:
        return None
    path = Path(path)
    if not path.is_absolute():
        path = Path(directory) / path
    return os.fspath(path)


def _spawn_detached_worker(directory, run_id):
    return subprocess.Popen(
        [
            sys.executable,
            "-m",
            "ase2sprkkr.execution.detached_worker",
            os.fspath(directory),
            run_id,
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
        close_fds=True,
    )


def _wait_for_direct_launch(worker, directory, run_id, poll_interval=0.01):
    poll = getattr(worker, "poll", None)
    if poll is None:
        return _read_state(directory, run_id)

    while True:
        state = _read_state(directory, run_id)
        if state.get("state") != RunState.SUBMITTING.value:
            return state

        return_code = poll()
        if return_code is None:
            time.sleep(poll_interval)
            continue

        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            if state.get("state") != RunState.SUBMITTING.value:
                return state
            return _write_state_unlocked(
                directory,
                run_id,
                state="failed",
                error=(
                    "Detached worker exited before calculation launch completed "
                    f"with return code {return_code}"
                ),
                failure_phase=FailurePhase.LAUNCH.value,
                return_code=None,
                finished_at=state.get("finished_at") or _timestamp(),
            )


def submit_detached(request):
    """Start one crash-surviving local calculation."""
    executor = request.get("executor", "local")
    if executor != "local":
        raise ValueError(
            f"Detached persistence supports only local execution, got {executor!r}"
        )

    directory = Path(request["directory"]).expanduser().resolve()
    durable_mkdir(directory, parents=True, exist_ok=True)
    run_id = request["id"]
    run_directory = _run_directory(directory, run_id)
    try:
        durable_mkdir(run_directory, parents=True)
    except FileExistsError as exc:
        raise RuntimeError(f"Detached run already exists: {run_id}") from exc

    request = dict(request, directory=os.fspath(directory), executor="local")
    input_file = _absolute_file(directory, request.get("input_file"))
    output_file = _absolute_file(directory, request.get("output_file"))
    error_file = _absolute_file(directory, request.get("error_file"))

    if input_file is not None and os.path.isfile(input_file):
        sync_file(input_file)
    durable_write_json(_request_path(directory, run_id), request)
    _write_state(
        directory,
        run_id,
        id=run_id,
        state=RunState.SUBMITTING.value,
        directory=os.fspath(directory),
        input_file=input_file,
        output_file=output_file,
        error_file=error_file,
        stream_output_file=None,
        stream_error_file=None,
        return_code=None,
        error=None,
        failure_phase=None,
        pid=None,
        created_at=_timestamp(),
        submitted_at=None,
        started_at=None,
        finished_at=None,
    )

    try:
        worker = _spawn_detached_worker(directory, run_id)
    except Exception as exc:
        return _write_state(
            directory,
            run_id,
            state="failed",
            error=repr(exc),
            failure_phase=FailurePhase.LAUNCH.value,
            return_code=None,
            finished_at=_timestamp(),
        )

    return _wait_for_direct_launch(worker, directory, run_id)


def status_detached(directory, run_id):
    return _read_state(Path(directory).expanduser().resolve(), run_id)


def _process_group_running(pid):
    try:
        os.killpg(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _wait_process_group_exit(pid, timeout):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not _process_group_running(pid):
            return True
        time.sleep(0.05)
    return not _process_group_running(pid)


def _terminate_process_group(directory, run_id, pid, timeout):
    if not pid:
        return _read_state(directory, run_id), True, False

    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return _read_state(directory, run_id), True, False

    if _wait_process_group_exit(pid, timeout):
        return _read_state(directory, run_id), True, True

    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        return _read_state(directory, run_id), True, True

    stopped = _wait_process_group_exit(pid, timeout)
    return _read_state(directory, run_id), stopped, True


def cancel_detached(directory, run_id, timeout=1.0):
    directory = Path(directory).expanduser().resolve()
    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATES:
            return state
        pid = state.get("pid")
        if not pid:
            return _write_state_unlocked(
                directory,
                run_id,
                state="cancelled",
                failure_phase=None,
                finished_at=_timestamp(),
            )

    state, stopped, signalled = _terminate_process_group(
        directory, run_id, pid, timeout
    )
    if not stopped or not signalled:
        return state

    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATES:
            return state
        return _write_state_unlocked(
            directory,
            run_id,
            state="cancelled",
            failure_phase=None,
            finished_at=state.get("finished_at") or _timestamp(),
        )


def follow_detached(
    directory,
    run_id,
    stdout,
    stderr,
    *,
    stdout_offset=0,
    stderr_offset=0,
    poll_interval=0.1,
):
    directory = Path(directory).expanduser().resolve()
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}
    while True:
        state = _read_state(directory, run_id)
        had_data = False
        seen_paths = set()
        for kind, key in (("out", "output_file"), ("err", "error_file")):
            path = detached_stream_file(state, key)
            if not path or path in seen_paths:
                continue
            seen_paths.add(path)
            try:
                with open(path, "rb") as source:
                    source.seek(offsets[kind])
                    data = source.read()
            except FileNotFoundError:
                data = b""
            if data:
                streams[kind].write(data)
                streams[kind].flush()
                offsets[kind] += len(data)
                had_data = True

        if state.get("state") in _FINAL_STATES and not had_data:
            return state
        if not had_data:
            time.sleep(poll_interval)


def list_detached_runs(directory):
    directory = Path(directory).expanduser().resolve()
    runs_directory = directory / _METADATA_DIR / _RUNS_DIR
    if not runs_directory.is_dir():
        return []

    runs = []
    for run_directory in runs_directory.iterdir():
        if not run_directory.is_dir():
            continue
        try:
            runs.append(_read_state(directory, run_directory.name))
        except (FileNotFoundError, json.JSONDecodeError):
            continue
    return sorted(runs, key=lambda state: state.get("created_at") or "")


def clean_detached_runs(directory, older_than=None):
    directory = Path(directory).expanduser().resolve()
    now = datetime.now(timezone.utc)
    removed = []
    active = []

    for state in list_detached_runs(directory):
        run_id = state["id"]
        if state.get("state") not in _FINAL_STATES:
            active.append(run_id)
            continue
        if older_than is not None:
            timestamp = state.get("finished_at") or state.get("created_at")
            if timestamp is None:
                continue
            finished = datetime.fromisoformat(timestamp)
            if finished.tzinfo is None:
                finished = finished.replace(tzinfo=timezone.utc)
            if (now - finished).total_seconds() < older_than:
                continue
        shutil.rmtree(_run_directory(directory, run_id))
        removed.append(run_id)

    return {"removed": removed, "active": active}


def resolve_detached_mpi(mpi):
    return LaunchSpec.from_mpi(mpi)


def _request_spec(directory, request):
    mpi_argv = request.get("mpi_argv")
    return ExecutionSpec(
        argv=tuple(request["argv"]),
        input_file=request.get("input_file"),
        output_file=request.get("output_file"),
        error_file=request.get("error_file"),
        directory=os.fspath(directory),
        stdin=request.get("stdin"),
        mpi_argv=None if mpi_argv is None else tuple(mpi_argv),
        launch=LaunchSpec(**(request.get("launch") or {})),
    )


def run_detached_worker(directory, run_id):
    directory = Path(directory).expanduser().resolve()
    try:
        request = json.loads(_request_path(directory, run_id).read_text())
        if request.get("executor", "local") != "local":
            raise ValueError("Detached worker received a non-local request")
        from .detached_local import run_detached_local

        run_detached_local(directory, run_id, request)
    except BaseException as exc:
        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            if state.get("state") != "cancelled":
                phase = (
                    FailurePhase.LAUNCH
                    if state.get("state") in {"submitting", "queued"}
                    else state.get("failure_phase")
                    or FailurePhase.EXECUTION.value
                )
                _write_state_unlocked(
                    directory,
                    run_id,
                    state="failed",
                    error=repr(exc),
                    failure_phase=(phase.value if isinstance(phase, FailurePhase) else phase),
                    return_code=state.get("return_code"),
                    finished_at=state.get("finished_at") or _timestamp(),
                )
        raise
