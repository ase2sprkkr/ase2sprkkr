"""Stateless target-side operations for remote execution.

The remote host is deliberately not a lifecycle authority. Every operation is
fully described by its request or by an opaque backend handle returned from
submission. No ASE2SPRKKR run record is written on the target.
"""

from __future__ import annotations

import asyncio
import base64
import json
import os
import re
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path
from uuid import uuid4

from .growing_file_tail import follow_growing_files
from .file_availability import clear_declared_outputs
from .base import ExecutionSpec, FailurePhase, LaunchSpec, RunState, SubmissionNotStarted
from .remote_process_worker import _SPEC_ENV, _START_FD_ENV
from .remote_mailbox import publish_outcome, read_outcome, acknowledge_outcome, _path as _outcome_path
from ._durability import durable_mkdir, durable_write_json, fsync_directory


_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _submission_ticket(run_id):
    if not isinstance(run_id, str) or re.fullmatch(r"[0-9a-f]{32}", run_id) is None:
        raise ValueError("Invalid project run id")
    return _outcome_path(run_id).parent.parent / "remote-submissions" / (run_id + ".json")


def _publish_submission(run_id, reply):
    """Temporary identity handoff, never a target-side run lifecycle record."""
    path = _submission_ticket(run_id)
    durable_mkdir(path.parent, parents=True, exist_ok=True, mode=0o700)
    durable_write_json(path, reply, mode=0o600)


def _recover_submission(run_id):
    try:
        return json.loads(_submission_ticket(run_id).read_text())
    except FileNotFoundError:
        return None


def acknowledge_submission(run_id):
    path = _submission_ticket(run_id)
    try:
        path.unlink()
    except FileNotFoundError:
        return {"acknowledged": False}
    fsync_directory(path.parent)
    return {"acknowledged": True}


def _encode_handle(data):
    payload = json.dumps(data, separators=(",", ":"), sort_keys=True).encode("utf8")
    return base64.urlsafe_b64encode(payload).decode("ascii").rstrip("=")


def _decode_handle(value):
    value = str(value)
    padding = "=" * (-len(value) % 4)
    try:
        data = json.loads(base64.urlsafe_b64decode(value + padding).decode("utf8"))
    except (ValueError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid remote execution id") from exc
    if not isinstance(data, dict) or "kind" not in data:
        raise ValueError("Invalid remote execution id")
    return data


def _absolute(directory, value):
    if value is None:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = Path(directory) / path
    return os.fspath(path)


def _process_info(pid):
    try:
        data = Path(f"/proc/{int(pid)}/stat").read_text()
    except (FileNotFoundError, ProcessLookupError, ValueError):
        return None
    end = data.rfind(")")
    if end < 0:
        return None
    fields = data[end + 2 :].split()
    if len(fields) <= 19:
        return None
    return fields[19], fields[0]


def _process_running(handle):
    info = _process_info(handle["pid"])
    return (
        info is not None
        and info[0] == str(handle["start_time"])
        and info[1] != "Z"
    )


def _request_spec(request):
    mpi_argv = request.get("mpi_argv")
    return ExecutionSpec(
        argv=tuple(request["argv"]),
        input_file=request.get("input_file"),
        output_file=request.get("output_file"),
        error_file=request.get("error_file"),
        directory=os.fspath(request["directory"]),
        stdin=request.get("stdin"),
        mpi_argv=None if mpi_argv is None else tuple(mpi_argv),
        launch=LaunchSpec(**(request.get("launch") or {})),
    )


def _local_submit(request):
    from .local import materialize_local_spec

    # Everything here happens before launching a supervisor. A failed local
    # filesystem setup proves that submission never started, unlike a lost
    # reply after Popen.
    try:
        directory = Path(request["directory"]).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        spec = materialize_local_spec(
            _request_spec(dict(request, directory=os.fspath(directory)))
        )

        output_path = _absolute(directory, spec.output_file)
        error_path = _absolute(directory, spec.error_file)
        for path in (output_path, error_path):
            if path is not None:
                Path(path).parent.mkdir(parents=True, exist_ok=True)
        sync = request.get("client_sync") or {}
        declared = (
            *request.get("output_files", ()),
            *sync.get("output_files", ()),
            *sync.get("output_patterns", ()),
            output_path, error_path,
        )
        protected = (
            request.get("input_file"), request.get("stdin"),
            *request.get("input_files", ()),
            *request.get("cache_files", ()),
            *(item.get("path") for item in sync.get("cache_files", ())),
            *sync.get("cache_patterns", ()),
        )
        clear_declared_outputs(declared, protected=protected)
    except OSError as exc:
        raise SubmissionNotStarted(str(exc)) from exc

    run_id = str(request.get("id") or uuid4().hex)
    receipt_id = uuid4().hex

    worker_spec = {
        "argv": list(spec.argv),
        "directory": os.fspath(directory),
        "stdin": _absolute(directory, spec.stdin),
        "output_file": output_path,
        "error_file": error_path,
        "receipt_id": receipt_id,
        "detached": bool(request.get("detached", True)),
    }
    env = os.environ.copy()
    env["ASE2SPRKKR_RUN_ID"] = run_id
    
    start_fd, release_fd = os.pipe()
    socket_read_fd, socket_write_fd = os.pipe()
    worker_spec["socket_fd"] = socket_write_fd
    env[_SPEC_ENV] = json.dumps(worker_spec, separators=(",", ":"))
    env[_START_FD_ENV] = str(start_fd)
    try:
        try:
            try:
                process = subprocess.Popen(
                    [sys.executable, "-m", "ase2sprkkr.execution.remote_process_worker"],
                    cwd=directory,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    env=env,
                    start_new_session=True,
                    close_fds=True,
                    pass_fds=(start_fd, socket_write_fd),
                )
            except OSError as exc:
                raise SubmissionNotStarted(str(exc)) from exc
        finally:
            os.close(start_fd)
            os.close(socket_write_fd)

        info = _process_info(process.pid)
        if info is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            raise RuntimeError("Started remote process disappeared before identification")
        start_time, _ = info
        os.write(release_fd, b"1")
        # An advertised handle must refer to an already-listening supervisor.
        with os.fdopen(socket_read_fd, "rb", closefd=True) as ready:
            socket_path = ready.readline().decode("utf8").rstrip("\n")
        if not socket_path:
            raise RuntimeError("Remote supervisor did not create its follow socket")
    finally:
        os.close(release_fd)

    remote_id = _encode_handle(
        {
            "kind": "process",
            "pid": process.pid,
            "start_time": start_time,
            "receipt_id": receipt_id,
            "socket_path": socket_path,
        }
    )
    reply = {
        "remote_id": remote_id,
        "state": RunState.RUNNING.value,
        "directory": os.fspath(directory),
        "input_file": _absolute(directory, spec.input_file),
        "output_file": output_path,
        "error_file": error_path,
        "return_code": None,
        "failure_phase": None,
    }
    if request.get("detached"):
        _publish_submission(run_id, {"remote_id": remote_id})
    return reply


def _psij_executor(options):
    from .psij import PsijExecutor

    options = options or {}
    return PsijExecutor(options.get("backend") or "slurm", url=options.get("url"))


def _psij_submit(request):
    executor = _psij_executor(request.get("executor_options"))
    spec = _request_spec(request)
    runner = (
        executor.create_project_run(spec, run_id=request["id"])
        if request.get("id")
        else executor.create_run(spec)
    )
    runner.submit()
    backend_id = runner.backend_id
    if backend_id is None:
        raise RuntimeError("PSI/J submission did not return a backend id")
    remote_id = _encode_handle(
        {
            "kind": "psij",
            "backend": executor.backend,
            "url": executor.url,
            "native_id": backend_id,
        }
    )
    return {
        "remote_id": remote_id,
        "state": runner.state.value,
        "directory": request["directory"],
        "input_file": request.get("input_file"),
        "output_file": runner.result.output_file,
        "error_file": runner.result.error_file,
        "return_code": runner.return_code,
        "failure_phase": (
            None if runner.failure_phase is None else runner.failure_phase.value
        ),
    }


def submit_remote(request):
    executor = request.get("executor", "local")
    if executor == "local":
        return _local_submit(request)
    if executor == "psij":
        return _psij_submit(request)
    raise ValueError(f"Unsupported remote executor: {executor!r}")


def _process_outcome(receipt_id):
    data = read_outcome(receipt_id)
    if data is None:
        return None
    try:
        state = RunState(data["state"])
        return_code = data.get("return_code")
        raw_phase = data.get("failure_phase")
        failure_phase = None if raw_phase is None else FailurePhase(raw_phase)
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("Invalid remote process outcome artifact") from exc
    return state, return_code, failure_phase, data.get("error")


def status_remote(remote_id):
    """Query backend facts; the client supplies any terminal receipt path."""
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        outcome = _process_outcome(handle["receipt_id"])
        if outcome is not None:
            state, return_code, failure_phase, error = outcome
        else:
            state = RunState.RUNNING if _process_running(handle) else RunState.UNKNOWN
            return_code = None
            failure_phase = None
            error = None
        return {
            "state": state.value,
            "return_code": return_code,
            "error": error,
            "failure_phase": (
                None if failure_phase is None else failure_phase.value
            ),
        }

    if handle["kind"] == "psij":
        if str(handle["backend"]).lower() == "slurm" and handle.get("url") is None:
            from .psij import _slurm_status_native
            status = _slurm_status_native(handle["native_id"])
        else:
            executor = _psij_executor(
                {"backend": handle["backend"], "url": handle.get("url")}
            )
            status = executor.status_native(handle["native_id"])
        return {
            "state": status.state.value,
            "return_code": status.return_code,
            "error": status.message,
            "failure_phase": (
                None if status.failure_phase is None else status.failure_phase.value
            ),
        }
    raise ValueError("Unsupported remote execution id")


def cancel_remote(remote_id, timeout=1.0):
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        outcome = _process_outcome(handle["receipt_id"])
        if outcome is not None:
            return status_remote(remote_id)
        try:
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as conn:
                conn.settimeout(timeout)
                conn.connect(handle["socket_path"])
                conn.sendall(b'{"command":"cancel"}\n')
                reply = conn.recv(256)
                if b'"accepted":true' not in reply:
                    raise RuntimeError("Supervisor rejected cancellation")
        except (FileNotFoundError, ConnectionRefusedError):
            return status_remote(remote_id)
        # The request is not a terminal fact. Only the supervisor waits for and
        # publishes the final child outcome.
        return status_remote(remote_id)

    if handle["kind"] == "psij":
        executor = _psij_executor(
            {"backend": handle["backend"], "url": handle.get("url")}
        )
        job = executor._psij.Job()
        executor._job_executor.attach(job, handle["native_id"])
        job.cancel()
        state = status_remote(remote_id)
        return state
    raise ValueError("Unsupported remote execution id")


def acknowledge_remote(remote_id):
    """Drop the terminal receipt only after the client persists its run state."""
    handle = _decode_handle(remote_id)
    if handle["kind"] != "process":
        return {"acknowledged": False}
    return {"acknowledged": acknowledge_outcome(handle["receipt_id"])}


def recover_remote(run_id, request, *, submitted_at=None):
    """Recover an uncertain submit using backend-native identity only."""
    executor = request.get("executor", "local")
    if executor == "local":
        reply = _recover_submission(run_id)
        if reply is None:
            return None
        return dict(reply, directory=request["directory"], **status_remote(reply["remote_id"]))
    if executor != "psij":
        raise ValueError(f"Unsupported remote executor: {executor!r}")

    backend = _psij_executor(request.get("executor_options"))
    facts = backend.recover_native(
        run_id, request["directory"], submitted_at=submitted_at,
    )
    if facts is None:
        return None
    native_id = facts["native_id"]
    remote_id = _encode_handle({
        "kind": "psij", "backend": backend.backend,
        "url": backend.url, "native_id": native_id,
    })
    state = facts["state"]
    return {
        "remote_id": remote_id,
        "state": state.value,
        "return_code": facts["return_code"],
        "failure_phase": (
            FailurePhase.EXECUTION.value if state == RunState.FAILED else None
        ),
        "directory": facts["directory"],
        "output_file": facts["stdout"],
        "error_file": facts["stderr"],
    }


async def _follow_scheduler_async(remote_id, streams, paths, offsets, poll_interval):
    """Follow scheduler files directly: no process supervisor is involved."""
    def emit(kind, data):
        streams[kind].write(data)
        streams[kind].flush()

    async def terminal():
        status = await asyncio.to_thread(status_remote, remote_id)
        return status if status["state"] in {s.value for s in _FINAL_STATES} else None

    return await follow_growing_files(
        paths, emit, terminal, offsets=offsets,
        poll_interval=max(poll_interval, 1.0),
    )


def _follow_scheduler(remote_id, streams, paths, offsets, poll_interval):
    return asyncio.run(
        _follow_scheduler_async(remote_id, streams, paths, offsets, poll_interval)
    )


def follow_remote(
    remote_id,
    stdout,
    stderr,
    *,
    stdout_offset=0,
    stderr_offset=0,
    poll_interval=0.1,
    stdout_path=None,
    stderr_path=None,
):
    """Follow direct supervisors or scheduler-owned streams until terminal state."""
    handle = _decode_handle(remote_id)
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}
    paths = {
        "out": stdout_path,
        "err": stderr_path,
    }

    if handle["kind"] == "psij":
        return _follow_scheduler(remote_id, streams, paths, offsets, poll_interval)
    if handle["kind"] != "process":
        raise ValueError("Unknown execution backend")
    from .remote_supervisor_stream import follow_socket
    path = handle.get("socket_path")
    if path is None:
        raise RuntimeError("Backend handle contains no follow socket")

    # Use the same frame decoder for the local socket and the remote CLI.
    # A merged stream has only one logical channel.
    def on_chunk(kind, offset, data):
        previous = offsets[kind]
        if previous is not None and offset != previous and paths[kind] is not None:
            raise RuntimeError("Supervisor FOLLOW offset mismatch")
        offsets[kind] = offset + len(data)
        streams[kind].write(data)
        streams[kind].flush()

    async def drain_from_socket():
        requested = dict(offsets)
        if stdout_path is not None and stdout_path == stderr_path:
            requested.pop("err")
        return await follow_socket(path, requested, on_chunk)

    try:
        return asyncio.run(drain_from_socket())
    except (FileNotFoundError, ConnectionRefusedError):
        # A completed supervisor has removed its transient socket. Replay
        # declared spools once, without the historical polling lifecycle.
        status = status_remote(remote_id)
        if status["state"] not in {s.value for s in _FINAL_STATES}:
            raise
        if not any(path is not None for path in paths.values()):
            raise ValueError("Output was live-only and cannot be reattached")
        seen = set()
        for kind in ("out", "err"):
            path = paths[kind]
            if path is None or path in seen:
                continue
            seen.add(path)
            with open(path, "rb") as source:
                source.seek(offsets[kind])
                while True:
                    data = source.read(65536)
                    if not data:
                        break
                    streams[kind].write(data)
                    streams[kind].flush()
        return status
