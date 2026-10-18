"""Stateless target-side operations for remote execution.

The remote host is deliberately not a lifecycle authority. Every operation is
fully described by its request or by an opaque backend handle returned from
submission. No ASE2SPRKKR run record is written on the target.
"""

from __future__ import annotations

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

from .base import ExecutionSpec, FailurePhase, LaunchSpec, RunState
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
    runner = executor.create_persistent_run(
        _request_spec(request),
        run_id=request.get("id"),
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


def _slurm_accounting_status(native_id):
    """Return exact scheduler facts for one Slurm id after it left scontrol."""
    from .psij import _scheduler_output, _slurm_exit_code, _slurm_state

    try:
        output = _scheduler_output(
            (
                "sacct",
                "-X",
                "-n",
                "-P",
                "-j",
                str(native_id),
                "--format=JobIDRaw,State,ExitCode",
            )
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None

    for line in output.splitlines():
        fields = line.rstrip("|").split("|")
        if len(fields) != 3 or fields[0].strip() != str(native_id):
            continue
        state = _slurm_state(fields[1])
        return state, _slurm_exit_code(fields[2].strip())
    return None


def _slurm_psij_status(handle):
    """Query Slurm directly; no PSI/J object may be target-side state memory."""
    from .psij import _slurm_details, _slurm_exit_code, _slurm_state

    native_id = str(handle["native_id"])
    try:
        details = _slurm_details(native_id)
    except (FileNotFoundError, subprocess.CalledProcessError):
        details = None

    if details:
        state = _slurm_state(details.get("JobState"))
        return_code = _slurm_exit_code(details.get("ExitCode"))
    else:
        accounted = _slurm_accounting_status(native_id)
        if accounted is None:
            state = RunState.UNKNOWN
            return_code = None
        else:
            state, return_code = accounted

    if state not in _FINAL_STATES:
        return_code = None
    elif state == RunState.FINISHED and return_code is None:
        return_code = 0
    failure_phase = FailurePhase.EXECUTION if state == RunState.FAILED else None
    return {
        "state": state.value,
        "return_code": return_code,
        "error": None,
        "failure_phase": (
            None if failure_phase is None else failure_phase.value
        ),
    }


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
        # PSI/J attach is asynchronous: a freshly attached Job initially reports
        # NEW and only later gets its real state. A stateless CLI process cannot
        # keep that Job object between requests. Slurm therefore has to be
        # queried by its native id instead of using PSI/J as target-side memory.
        if str(handle.get("backend", "")).lower() == "slurm":
            return _slurm_psij_status(handle)
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
        if state["state"] not in {item.value for item in _FINAL_STATES}:
            state["state"] = RunState.CANCELLED.value
            state["failure_phase"] = None
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


def _follow_scheduler(remote_id, streams, paths, offsets, poll_interval):
    """Scheduler-native stream observation; no direct-process fallback.

    A PSI/J scheduler does not run the direct-process supervisor. Its output
    is captured into scheduler-managed files and its terminal status comes
    from the native backend. Read at most one bounded chunk per channel.
    """
    while True:
        had_data = False
        seen = set()
        for kind in ("out", "err"):
            path = paths[kind]
            if path is None or path in seen:
                continue
            seen.add(path)
            try:
                with open(path, "rb") as source:
                    if offsets[kind] is None:
                        source.seek(0, os.SEEK_END)
                        offsets[kind] = source.tell()
                    else:
                        source.seek(offsets[kind])
                    data = source.read(65536)
            except FileNotFoundError:
                data = b""
            if data:
                streams[kind].write(data)
                streams[kind].flush()
                offsets[kind] += len(data)
                had_data = True
        status = status_remote(remote_id)
        if status["state"] in {state.value for state in _FINAL_STATES} | {RunState.UNKNOWN.value}:
            if not had_data:
                return status
        if not had_data:
            time.sleep(poll_interval)


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
    from .remote_supervisor_stream import MAGIC, _HEADER
    path = handle.get("socket_path")
    if path is None:
        raise RuntimeError("Backend handle contains no follow socket")

    def drain_from_socket():
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
            sock.connect(path)
            offsets = {"out": stdout_offset, "err": stderr_offset}
            if stdout_path is not None and stdout_path == stderr_path:
                offsets.pop("err")
            sock.sendall(json.dumps({"command": "follow", "offsets": offsets}).encode() + b"\n")

            def read_exactly(n):
                parts = []
                while n:
                    part = sock.recv(n)
                    if not part:
                        raise ConnectionError("Supervisor FOLLOW ended before terminal state")
                    parts.append(part)
                    n -= len(part)
                return b"".join(parts)

            if read_exactly(len(MAGIC)) != MAGIC:
                raise RuntimeError("Unknown supervisor protocol")
            while True:
                channel, offset, size = _HEADER.unpack(read_exactly(_HEADER.size))
                if channel == b"T":
                    return json.loads(read_exactly(size))
                if channel not in (b"O", b"E") or size > 65536:
                    raise RuntimeError("Invalid supervisor FOLLOW frame")
                kind = "out" if channel == b"O" else "err"
                data = read_exactly(size)
                # Live-only streams are intentionally non-replayable and may
                # start after bytes have already been emitted.
                if offsets[kind] is not None and offset != offsets[kind] and paths[kind] is not None:
                    raise RuntimeError("Supervisor FOLLOW offset mismatch")
                offsets[kind] = offset + len(data)
                streams[kind].write(data)
                streams[kind].flush()

    try:
        return drain_from_socket()
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
