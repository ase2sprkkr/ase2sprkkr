"""Shared detached execution protocol for local and remote persistent runs.

This module owns the target-side state files and worker lifecycle.  It is
independent of how a target is reached: local persistent execution uses it on
the current machine, while remote execution invokes the same protocol through
a transport.
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
from .base import ExecutionSpec, LaunchSpec


_FINAL_STATES = {"finished", "failed", "cancelled"}
_METADATA_DIR = ".ase2sprkkr"
_RUNS_DIR = "runs"
_STATE_FILE = "state.json"
_STATE_LOCK_FILE = "state.lock"
_REQUEST_FILE = "request.json"
_DETACHED_EXECUTORS = {"local", "psij"}


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
    """Return the inter-process lock guarding one detached run state."""
    return FileLock(_state_lock_path(directory, run_id))


def _read_state_unlocked(directory, run_id):
    with _state_path(directory, run_id).open() as file:
        return json.load(file)


def _read_state(directory, run_id):
    return _read_state_unlocked(directory, run_id)


def _write_state_unlocked(workspace, run_id, **changes):
    path = _state_path(workspace, run_id)
    try:
        state = _read_state_unlocked(workspace, run_id)
    except FileNotFoundError:
        state = {}
    state.update(changes)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(state))
    os.replace(temporary, path)
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


def submit_detached(request):
    """Start a detached target-side executor and return its initial state."""
    executor = request.get("executor", "local")
    if executor not in _DETACHED_EXECUTORS:
        raise ValueError(f"Unsupported detached executor: {executor!r}")

    directory = Path(request["directory"]).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    run_id = request["id"]
    run_directory = _run_directory(directory, run_id)
    try:
        run_directory.mkdir(parents=True)
    except FileExistsError as exc:
        raise RuntimeError(f"Detached run already exists: {run_id}") from exc

    request = dict(request, directory=os.fspath(directory))
    output_file = _absolute_file(directory, request.get("output_file"))
    error_file = _absolute_file(directory, request.get("error_file"))
    cleared = set()
    for path in (output_file, error_file):
        if path is None or path in cleared:
            continue
        cleared.add(path)
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"")

    _request_path(directory, run_id).write_text(json.dumps(request))
    state = _write_state(
        directory,
        run_id,
        id=run_id,
        state="queued",
        directory=os.fspath(directory),
        input_file=_absolute_file(directory, request.get("input_file")),
        output_file=output_file,
        error_file=error_file,
        return_code=None,
        error=None,
        pid=None,
        backend_reference=None,
        created_at=_timestamp(),
        started_at=None,
        finished_at=None,
    )

    subprocess.Popen(
        [
            sys.executable,
            "-m",
            "ase2sprkkr.execution.remote_worker",
            os.fspath(directory),
            run_id,
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
        close_fds=True,
    )
    return state


def status_detached(directory, run_id):
    return _read_state(Path(directory).expanduser().resolve(), run_id)


def _terminate_process_group(directory, run_id, pid, timeout):
    if not pid:
        return _read_state(directory, run_id)

    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return _read_state(directory, run_id)

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.killpg(pid, 0)
        except ProcessLookupError:
            return _read_state(directory, run_id)
        except PermissionError:
            pass
        time.sleep(0.05)

    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return _read_state(directory, run_id)


def cancel_detached(directory, run_id, timeout=1.0):
    """Cancel a detached run through its target-side executor."""
    directory = Path(directory).expanduser().resolve()
    request = json.loads(_request_path(directory, run_id).read_text())
    if request.get("executor", "local") == "psij":
        from .remote_psij import cancel_remote_psij

        return cancel_remote_psij(directory, run_id, request, timeout=timeout)

    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATES:
            return state
        state = _write_state_unlocked(
            directory, run_id, state="cancelled", finished_at=_timestamp()
        )
        pid = state.get("pid")

    return _terminate_process_group(directory, run_id, pid, timeout)


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
    """Forward newly written stdout/stderr until the detached run is final.

    If both state entries name the same physical file, it is an explicitly
    merged stream. Forward it once as stdout; the original stdout/stderr
    distinction no longer exists after such a merge.
    """
    directory = Path(directory).expanduser().resolve()
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}

    while True:
        state = _read_state(directory, run_id)
        had_data = False
        seen_paths = set()
        for kind, key in (("out", "output_file"), ("err", "error_file")):
            path = state.get(key)
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
    """Return states of detached runs recorded in a workspace."""
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
    """Remove completed detached metadata without touching workspace artifacts."""
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
    """Compatibility helper for the detached target-side CLI."""
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
    """Execute one run using the executor selected by the detached request."""
    from .local import materialize_local_spec

    directory = Path(directory).expanduser().resolve()
    try:
        request = json.loads(_request_path(directory, run_id).read_text())
        executor = request.get("executor", "local")
        if executor == "psij":
            from .remote_psij import run_remote_psij_worker

            run_remote_psij_worker(directory, run_id, request)
            return
        if executor != "local":
            raise ValueError(f"Unsupported detached executor: {executor!r}")

        spec = materialize_local_spec(_request_spec(directory, request))
        output_path = Path(spec.output_file or "stdout")
        error_path = Path(spec.error_file or "stderr")
        if not output_path.is_absolute():
            output_path = directory / output_path
        if not error_path.is_absolute():
            error_path = directory / error_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        error_path.parent.mkdir(parents=True, exist_ok=True)

        if spec.stdin:
            stdin_path = Path(spec.stdin)
            if not stdin_path.is_absolute():
                stdin_path = directory / stdin_path
            stdin = stdin_path.open("rb")
        else:
            stdin = subprocess.DEVNULL

        def run_process(output, error):
            with _state_lock(directory, run_id):
                state = _read_state_unlocked(directory, run_id)
                if state.get("state") in _FINAL_STATES:
                    return None
                if state.get("state") != "queued":
                    raise RuntimeError(
                        f"Cannot start detached run {run_id!r} from state "
                        f"{state.get('state')!r}"
                    )
                process = subprocess.Popen(
                    spec.argv,
                    cwd=directory,
                    stdin=stdin,
                    stdout=output,
                    stderr=error,
                    start_new_session=True,
                )
                _write_state_unlocked(
                    directory,
                    run_id,
                    state="running",
                    pid=process.pid,
                    started_at=_timestamp(),
                )
            return process.wait()

        try:
            with output_path.open("wb") as output:
                if error_path == output_path:
                    return_code = run_process(output, subprocess.STDOUT)
                else:
                    with error_path.open("wb") as error:
                        return_code = run_process(output, error)
        finally:
            if stdin is not subprocess.DEVNULL:
                stdin.close()

        if return_code is None:
            return

        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            final_state = state.get("state")
            if final_state != "cancelled":
                final_state = "finished" if return_code == 0 else "failed"
            _write_state_unlocked(
                directory,
                run_id,
                state=final_state,
                return_code=return_code,
                finished_at=state.get("finished_at") or _timestamp(),
            )
    except BaseException as exc:
        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            if state.get("state") != "cancelled":
                _write_state_unlocked(
                    directory,
                    run_id,
                    state="failed",
                    error=repr(exc),
                    return_code=None,
                    finished_at=_timestamp(),
                )
        raise
