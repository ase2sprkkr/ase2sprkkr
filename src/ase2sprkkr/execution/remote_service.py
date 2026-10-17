"""Stateless target-side operations for remote execution.

The remote host is deliberately not a lifecycle authority. Every operation is
fully described by its request or by an opaque backend handle returned from
submission. No ASE2SPRKKR run record is written on the target.
"""

from __future__ import annotations

import base64
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

from .base import ExecutionSpec, FailurePhase, LaunchSpec, RunReference, RunState
from .remote_process_worker import _SPEC_ENV, _START_FD_ENV, _write_outcome


_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


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


def _same_process(handle):
    info = _process_info(handle["pid"])
    return info is not None and info[0] == str(handle["start_time"])


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

    output_path = _absolute(directory, spec.output_file or "stdout")
    error_path = _absolute(directory, spec.error_file or "stderr")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(error_path).parent.mkdir(parents=True, exist_ok=True)

    run_id = str(request.get("id") or uuid4().hex)
    outcome_file = _absolute(
        directory, request.get("outcome_file") or f".{run_id}.outcome.json"
    )
    Path(outcome_file).parent.mkdir(parents=True, exist_ok=True)
    try:
        os.unlink(outcome_file)
    except FileNotFoundError:
        pass

    worker_spec = {
        "argv": list(spec.argv),
        "directory": os.fspath(directory),
        "stdin": _absolute(directory, spec.stdin),
        "output_file": output_path,
        "error_file": error_path,
        "outcome_file": outcome_file,
    }
    env = os.environ.copy()
    env["ASE2SPRKKR_RUN_ID"] = run_id
    env[_SPEC_ENV] = json.dumps(worker_spec, separators=(",", ":"))

    start_fd, release_fd = os.pipe()
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
                pass_fds=(start_fd,),
            )
        finally:
            os.close(start_fd)

        info = _process_info(process.pid)
        if info is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            raise RuntimeError("Started remote process disappeared before identification")
        start_time, _ = info
        os.write(release_fd, b"1")
    finally:
        os.close(release_fd)

    remote_id = _encode_handle(
        {
            "kind": "process",
            "pid": process.pid,
            "start_time": start_time,
        }
    )
    return {
        "remote_id": remote_id,
        "state": RunState.RUNNING.value,
        "directory": os.fspath(directory),
        "input_file": _absolute(directory, spec.input_file),
        "output_file": output_path,
        "error_file": error_path,
        "return_code": None,
        "failure_phase": None,
    }


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
    reference = runner.reference
    if reference is None:
        raise RuntimeError("PSI/J submission did not return a backend id")
    remote_id = _encode_handle(
        {
            "kind": "psij",
            "backend": executor.backend,
            "url": executor.url,
            "native_id": reference.native_id,
        }
    )
    return {
        "remote_id": remote_id,
        "state": runner.state.value,
        "directory": reference.directory,
        "input_file": reference.input_file,
        "output_file": reference.output_file,
        "error_file": reference.error_file,
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


def _reference(handle):
    return RunReference(
        executor=f"psij:{handle['backend']}",
        native_id=str(handle["native_id"]),
        computer=handle.get("url"),
    )


def _process_outcome(path):
    try:
        data = json.loads(Path(path).read_text())
    except FileNotFoundError:
        return None
    except (KeyError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError("Invalid remote process outcome artifact") from exc
    try:
        state = RunState(data["state"])
        return_code = data.get("return_code")
        raw_phase = data.get("failure_phase")
        failure_phase = None if raw_phase is None else FailurePhase(raw_phase)
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("Invalid remote process outcome artifact") from exc
    return state, return_code, failure_phase


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


def status_remote(remote_id, *, outcome_file=None):
    """Query backend facts; the client supplies any terminal receipt path."""
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        if _process_running(handle):
            state = RunState.RUNNING
            return_code = None
            failure_phase = None
        else:
            outcome = None if outcome_file is None else _process_outcome(outcome_file)
            if outcome is None:
                state = RunState.UNKNOWN
                return_code = None
                failure_phase = None
            else:
                state, return_code, failure_phase = outcome
        return {
            "state": state.value,
            "return_code": return_code,
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
        status = executor.status(_reference(handle))[0]
        return {
            "state": status.state.value,
            "return_code": status.return_code,
            "error": status.message,
            "failure_phase": (
                None if status.failure_phase is None else status.failure_phase.value
            ),
        }
    raise ValueError("Unsupported remote execution id")


def cancel_remote(remote_id, timeout=1.0, *, outcome_file=None):
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        if outcome_file is None:
            raise ValueError("Direct-process cancellation requires its outcome path")
        if not _same_process(handle):
            return status_remote(remote_id, outcome_file=outcome_file)

        try:
            os.killpg(int(handle["pid"]), signal.SIGTERM)
        except ProcessLookupError:
            return status_remote(remote_id, outcome_file=outcome_file)

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline and _process_running(handle):
            time.sleep(0.05)
        if _process_running(handle):
            try:
                os.killpg(int(handle["pid"]), signal.SIGKILL)
            except ProcessLookupError:
                pass
            # SIGKILL is a request, not proof of process termination.
            deadline = time.monotonic() + timeout
            while time.monotonic() < deadline and _process_running(handle):
                time.sleep(0.05)
            if _process_running(handle):
                return status_remote(remote_id)

        # A natural completion racing with cancellation wins. Otherwise publish
        # cancellation as the durable backend outcome so future status queries
        # cannot regress to FAILED/UNKNOWN after the supervisor disappears.
        current = status_remote(remote_id, outcome_file=outcome_file)
        if current["state"] in (RunState.FINISHED.value, RunState.FAILED.value):
            return current
        _write_outcome(
            outcome_file,
            RunState.CANCELLED,
            current.get("return_code"),
            None,
        )
        return status_remote(remote_id, outcome_file=outcome_file)

    if handle["kind"] == "psij":
        executor = _psij_executor(
            {"backend": handle["backend"], "url": handle.get("url")}
        )
        job = executor._psij.Job()
        executor._job_executor.attach(job, handle["native_id"])
        job.cancel()
        state = status_remote(remote_id, outcome_file=outcome_file)
        if state["state"] not in {item.value for item in _FINAL_STATES}:
            state["state"] = RunState.CANCELLED.value
            state["failure_phase"] = None
        return state
    raise ValueError("Unsupported remote execution id")


def recover_remote(run_id, request, *, submitted_at=None):
    """Recover an uncertain submit using backend-native identity only."""
    executor = request.get("executor", "local")
    if executor == "local":
        raise NotImplementedError(
            "A direct remote process cannot recover a lost submit response from "
            "PID identity alone"
        )
    if executor != "psij":
        raise ValueError(f"Unsupported remote executor: {executor!r}")

    backend = _psij_executor(request.get("executor_options"))
    status = backend.recover_submission(
        run_id,
        SimpleNamespace(directory=request["directory"]),
        submitted_at=submitted_at,
    )
    if status is None:
        return None
    reference = status.reference
    remote_id = _encode_handle(
        {
            "kind": "psij",
            "backend": backend.backend,
            "url": backend.url,
            "native_id": reference.native_id,
        }
    )
    return {
        "remote_id": remote_id,
        "state": status.state.value,
        "return_code": status.return_code,
        "error": status.message,
        "failure_phase": (
            None if status.failure_phase is None else status.failure_phase.value
        ),
        "directory": reference.directory,
        "input_file": reference.input_file,
        "output_file": reference.output_file,
        "error_file": reference.error_file,
    }


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
    outcome_file=None,
):
    """Tail client-selected stream files until backend termination."""
    handle = _decode_handle(remote_id)
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}
    paths = {
        "out": stdout_path,
        "err": stderr_path,
    }

    while True:
        had_data = False
        seen = set()
        for kind in ("out", "err"):
            path = paths[kind]
            if not path or path in seen:
                continue
            seen.add(path)
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

        state = status_remote(remote_id, outcome_file=outcome_file)
        terminal = state["state"] in {item.value for item in _FINAL_STATES}
        unresolved = state["state"] == RunState.UNKNOWN.value
        if (terminal or unresolved) and not had_data:
            return state
        if not had_data:
            time.sleep(poll_interval)
