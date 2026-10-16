"""Stateless target-side operations for remote execution.

The remote host is deliberately not a lifecycle authority. Each operation is
fully described by its request or by an opaque ``remote_id`` returned from
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

from .base import ExecutionSpec, LaunchSpec, RunReference, RunState
from .remote_process_worker import _SPEC_ENV


_FINAL_STATES = {
    RunState.FINISHED,
    RunState.EXECUTION_FAILED,
    RunState.CANCELLED,
}


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
    """Return Linux process start time and state, or ``None`` when absent."""
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
    return_code_file = os.fspath(directory / f".{run_id}.returncode")
    try:
        os.unlink(return_code_file)
    except FileNotFoundError:
        pass

    worker_spec = {
        "argv": list(spec.argv),
        "directory": os.fspath(directory),
        "stdin": _absolute(directory, spec.stdin),
        "output_file": output_path,
        "error_file": error_path,
        "return_code_file": return_code_file,
    }
    env = os.environ.copy()
    env["ASE2SPRKKR_RUN_ID"] = run_id
    env[_SPEC_ENV] = json.dumps(worker_spec, separators=(",", ":"))
    process = subprocess.Popen(
        [sys.executable, "-m", "ase2sprkkr.execution.remote_process_worker"],
        cwd=directory,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=env,
        start_new_session=True,
        close_fds=True,
    )
    info = _process_info(process.pid)
    if info is None:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        raise RuntimeError("Started remote process disappeared before identification")
    start_time, _ = info

    remote_id = _encode_handle(
        {
            "kind": "process",
            "pid": process.pid,
            "start_time": start_time,
            "return_code_file": return_code_file,
            "directory": os.fspath(directory),
            "input_file": _absolute(directory, spec.input_file),
            "output_file": output_path,
            "error_file": error_path,
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
    }


class _DrainResult:
    return_code = None


class _DrainReader:
    print_output = False

    def create_result(self, output_file, input_file):
        return _DrainResult()

    async def parse_streams(self, stdout, stderr, result):
        return await stdout.read(), await stderr.read()

    def complete_result(self, result, output, error, return_code):
        result.return_code = return_code
        return result


def _psij_executor(options):
    from .psij import PsijExecutor

    options = options or {}
    return PsijExecutor(options.get("backend") or "slurm", url=options.get("url"))


def _psij_submit(request):
    executor = _psij_executor(request.get("executor_options"))
    runner = executor.create_persistent_run(
        _DrainReader(),
        _request_spec(request),
        run_id=request.get("id"),
    )
    # PSI/J owns scheduler submission. The helper process exits immediately
    # after the native reference is available; its daemon monitor thread is not
    # a persistent lifecycle authority and cannot outlive this command.
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
            "directory": reference.directory,
            "input_file": reference.input_file,
            "output_file": reference.output_file,
            "error_file": reference.error_file,
        }
    )
    return {
        "remote_id": remote_id,
        "state": RunState.QUEUED.value,
        "directory": reference.directory,
        "input_file": reference.input_file,
        "output_file": reference.output_file,
        "error_file": reference.error_file,
        "return_code": None,
    }


def submit_remote(request):
    """Submit one run without creating target-side ASE2SPRKKR lifecycle state."""
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
        directory=handle.get("directory"),
        input_file=handle.get("input_file"),
        output_file=handle.get("output_file"),
        error_file=handle.get("error_file"),
    )


def _process_return_code(handle):
    try:
        return int(Path(handle["return_code_file"]).read_text().strip())
    except FileNotFoundError:
        return None
    except (KeyError, ValueError) as exc:
        raise RuntimeError("Invalid remote process return-code artifact") from exc


def status_remote(remote_id):
    """Query backend facts for an opaque target id."""
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        if _process_running(handle):
            state = RunState.RUNNING
            return_code = None
        else:
            return_code = _process_return_code(handle)
            state = (
                RunState.UNKNOWN
                if return_code is None
                else RunState.FINISHED
                if return_code == 0
                else RunState.EXECUTION_FAILED
            )
        return {
            "state": state.value,
            "return_code": return_code,
            "directory": handle.get("directory"),
            "input_file": handle.get("input_file"),
            "output_file": handle.get("output_file"),
            "error_file": handle.get("error_file"),
        }
    if handle["kind"] == "psij":
        executor = _psij_executor(
            {"backend": handle["backend"], "url": handle.get("url")}
        )
        status = executor.status(_reference(handle))[0]
        return {
            "state": status.state.value,
            "return_code": status.return_code,
            "error": status.message,
            "directory": handle.get("directory"),
            "input_file": handle.get("input_file"),
            "output_file": handle.get("output_file"),
            "error_file": handle.get("error_file"),
        }
    raise ValueError("Unsupported remote execution id")


def cancel_remote(remote_id, timeout=1.0):
    handle = _decode_handle(remote_id)
    if handle["kind"] == "process":
        if _same_process(handle):
            try:
                os.killpg(int(handle["pid"]), signal.SIGTERM)
            except ProcessLookupError:
                pass
            deadline = time.monotonic() + timeout
            while time.monotonic() < deadline and _process_running(handle):
                time.sleep(0.05)
            if _process_running(handle):
                try:
                    os.killpg(int(handle["pid"]), signal.SIGKILL)
                except ProcessLookupError:
                    pass
        state = status_remote(remote_id)
        state["state"] = RunState.CANCELLED.value
        return state
    if handle["kind"] == "psij":
        executor = _psij_executor(
            {"backend": handle["backend"], "url": handle.get("url")}
        )
        runner = executor.attach(_reference(handle), _DrainReader(), start=False)
        runner.cancel()
        state = status_remote(remote_id)
        if state["state"] not in {item.value for item in _FINAL_STATES}:
            state["state"] = RunState.CANCELLED.value
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
            "directory": reference.directory,
            "input_file": reference.input_file,
            "output_file": reference.output_file,
            "error_file": reference.error_file,
        }
    )
    return {
        "remote_id": remote_id,
        "state": status.state.value,
        "return_code": status.return_code,
        "error": status.message,
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
):
    """Tail backend output until it is no longer active."""
    handle = _decode_handle(remote_id)
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}
    paths = {
        "out": handle.get("output_file"),
        "err": handle.get("error_file"),
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

        state = status_remote(remote_id)
        terminal = state["state"] in {item.value for item in _FINAL_STATES}
        lost_process = (
            handle["kind"] == "process"
            and state["state"] == RunState.UNKNOWN.value
        )
        if (terminal or lost_process) and not had_data:
            return state
        if not had_data:
            time.sleep(poll_interval)
