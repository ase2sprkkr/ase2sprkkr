"""Stateless local process backend for reconnectable execution.

The worker owns a calculation process and durable output publication.  The two
small immutable receipts are backend facts (process identity and terminal
outcome), not a second RunRecord or a lifecycle state machine.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from ._durability import durable_mkdir, durable_replace, durable_write_json, sync_stream
from .base import FailurePhase, RunState, RunStatus, SubmissionNotStarted


def _run_dir(directory, run_id):
    return Path(directory).resolve() / ".ase2sprkkr" / "runs" / run_id


def _boot_id():
    return Path("/proc/sys/kernel/random/boot_id").read_text().strip()


def _process_start(pid):
    """Linux kernel start ticks; avoid treating a reused PID as our worker."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        # The command name is parenthesized and may contain spaces or ')'.
        fields = stat[stat.rfind(")") + 2:].split()
        return fields[19], fields[0]  # starttime (field 22), state (field 3)
    except (OSError, IndexError, ValueError):
        return None


def _identity(run_id, pid=None):
    pid = os.getpid() if pid is None else pid
    found = _process_start(pid)
    if found is None:
        raise RuntimeError("Reconnectable local processes require Linux /proc")
    return {"run_id": run_id, "pid": pid, "start": found[0], "boot": _boot_id()}


def _backend_id(identity):
    return "local:" + ":".join(
        (identity["run_id"], str(identity["pid"]), identity["start"], identity["boot"])
    )


def _parse_id(identifier):
    try:
        prefix, run_id, pid, start, boot = identifier.split(":", 4)
        if prefix != "local" or not run_id or int(pid) <= 0 or not start or not boot:
            return None
        return {"run_id": run_id, "pid": int(pid), "start": start, "boot": boot}
    except (TypeError, ValueError, AttributeError):
        return None


def _receipt(directory, run_id, name):
    path = _run_dir(directory, run_id) / name
    try:
        with path.open(encoding="utf8") as source:
            return json.load(source)
    except FileNotFoundError:
        return None


def _alive(identity):
    try:
        if _boot_id() != identity["boot"]:
            return False
        found = _process_start(identity["pid"])
        return found is not None and found[0] == identity["start"] and found[1] != "Z"
    except OSError:
        return False


def query(directory, identifier):
    identity = _parse_id(identifier)
    if identity is None:
        return RunStatus(state=RunState.UNKNOWN, backend_id=identifier)
    outcome = _receipt(directory, identity["run_id"], "outcome.json")
    if outcome is not None:
        return RunStatus(
            state=(RunState.CANCELLED if outcome["cancelled"] else
                   RunState.FAILED if outcome["phase"] is not None or outcome["return_code"] != 0
                   else RunState.FINISHED),
            backend_id=identifier,
            return_code=outcome["return_code"],
            message=outcome.get("message"),
            failure_phase=outcome.get("phase"),
        )
    alive = _alive(identity)
    return RunStatus(
        state=RunState.RUNNING if alive else RunState.UNKNOWN,
        backend_id=identifier,
        message=None if alive else "Local worker disappeared without a terminal outcome",
    )


def recover(directory, run_id):
    identity = _receipt(directory, run_id, "identity.json")
    if identity is None:
        return None
    if identity["run_id"] != run_id:
        raise ValueError("Local backend identity does not match the run")
    return query(directory, _backend_id(identity))


def _absolute(directory, value):
    if value is None:
        return None
    path = Path(value)
    return str(path if path.is_absolute() else Path(directory) / path)


def stream_paths(directory, run_id, output, error, *, terminal=False):
    meta = _run_dir(directory, run_id)
    out = Path(output) if output else meta / "stdout"
    err = Path(error) if error else meta / "stderr"
    def launching(path, kind):
        return path.with_name(f".{path.name}.{run_id}.{kind}.launching")
    if terminal:
        return str(out), str(err)
    return str(launching(out, "out")), str(launching(err, "err")) if err != out else str(launching(out, "out"))


def observed_paths(directory, identifier, output, error):
    identity = _parse_id(identifier)
    if identity is None:
        return (None, None)
    outcome = _receipt(directory, identity["run_id"], "outcome.json")
    launch = stream_paths(directory, identity["run_id"], output, error)
    canonical = stream_paths(directory, identity["run_id"], output, error, terminal=True)
    # Never replay a stale canonical file from an earlier run after failed publication.
    published = (False, False) if outcome is None else outcome.get("published", (False, False))
    return tuple(canonical[i] if published[i] else launch[i] for i in (0, 1))


def _terminate(identity):
    if not _alive(identity):
        return
    # pidfd pins the kernel PID and eliminates reuse between check and signal.
    if hasattr(os, "pidfd_open") and hasattr(signal, "pidfd_send_signal"):
        try:
            fd = os.pidfd_open(identity["pid"])
            try:
                if _alive(identity):
                    signal.pidfd_send_signal(fd, signal.SIGTERM)
            finally:
                os.close(fd)
        except ProcessLookupError:
            pass
    elif _alive(identity):
        os.kill(identity["pid"], signal.SIGTERM)


def cancel(directory, identifier):
    identity = _parse_id(identifier)
    if identity is None:
        raise ValueError("Unrecognized local backend identity")
    if _receipt(directory, identity["run_id"], "outcome.json") is None:
        _terminate(identity)  # non-blocking; worker publishes a terminal receipt


def submit(directory, run_id, request):
    """Launch worker with immutable argv; no target-side request/state database."""
    directory = str(Path(directory).resolve())
    durable_mkdir(_run_dir(directory, run_id), parents=True)
    command = [
        sys.executable, "-m", "ase2sprkkr.execution.local_backend",
        directory, run_id, json.dumps(request),
    ]
    try:
        child = subprocess.Popen(
            command, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True, close_fds=True,
        )
    except OSError as exc:
        raise SubmissionNotStarted(str(exc)) from exc

    # The worker persists its identity *before* starting the calculation.
    # A lost reply can recover by run_id, without a second lifecycle store.
    while True:
        identity = _receipt(directory, run_id, "identity.json")
        if identity is not None:
            if identity["pid"] != child.pid or identity["run_id"] != run_id:
                raise RuntimeError("Mismatched local worker identity")
            return _backend_id(identity)
        code = child.poll()
        if code is not None:
            raise SubmissionNotStarted(
                f"Local worker exited before the durable handoff (exit {code})"
            )
        time.sleep(0.01)


def _worker(directory, run_id, request):
    meta = _run_dir(directory, run_id)
    identity = _identity(run_id)
    cancelled = False
    process = None
    def on_cancel(_signum, _frame):
        nonlocal cancelled
        cancelled = True
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    signal.signal(signal.SIGTERM, on_cancel)
    durable_write_json(meta / "identity.json", identity)

    code = None
    phase = None
    message = None
    out_path, err_path = stream_paths(
        directory, run_id, request.get("output_file"), request.get("error_file"),
    )
    out_final, err_final = stream_paths(
        directory, run_id, request.get("output_file"), request.get("error_file"), terminal=True,
    )
    merged = out_path == err_path
    started = False
    publishing = False
    published = [False, False]
    try:
        for path in {out_path, err_path, out_final, err_final}:
            durable_mkdir(Path(path).parent, parents=True, exist_ok=True)
        if not cancelled:
            from .base import ExecutionSpec, LaunchSpec
            from .local_process import materialize_local_spec
            spec = materialize_local_spec(ExecutionSpec(
                argv=request["argv"], directory=directory,
                mpi_argv=request.get("mpi_argv"), stdin=request.get("stdin"),
                launch=LaunchSpec(**request["launch"]),
            ))
            with open(out_path, "wb") as stdout:
                stderr = stdout if merged else open(err_path, "wb")
                try:
                    stdin = open(spec.stdin, "rb") if spec.stdin else subprocess.DEVNULL
                    try:
                        if not cancelled:
                            process = subprocess.Popen(
                                spec.argv, cwd=directory, stdin=stdin,
                                stdout=stdout, stderr=stderr, start_new_session=True,
                            )
                            started = True
                            if cancelled:
                                try:
                                    os.killpg(process.pid, signal.SIGTERM)
                                except ProcessLookupError:
                                    pass
                            deadline = None
                            while True:
                                if cancelled and deadline is None:
                                    deadline = time.monotonic() + 1
                                if cancelled and deadline is not None and time.monotonic() >= deadline:
                                    if process.poll() is None:
                                        try:
                                            os.killpg(process.pid, signal.SIGKILL)
                                        except ProcessLookupError:
                                            pass
                                try:
                                    code = process.wait(timeout=0.1)
                                    break
                                except subprocess.TimeoutExpired:
                                    continue
                    finally:
                        if hasattr(stdin, "close"):
                            stdin.close()
                    publishing = True
                    sync_stream(stdout)
                    if not merged:
                        sync_stream(stderr)
                finally:
                    if not merged:
                        stderr.close()
        if not cancelled:
            publishing = True
            if started:
                durable_replace(out_path, out_final)
                published[0] = True
                if not merged:
                    durable_replace(err_path, err_final)
                published[1] = True
    except BaseException as exc:
        phase = (FailurePhase.ARTIFACTS if publishing else
                 FailurePhase.EXECUTION if started else FailurePhase.LAUNCH)
        message = repr(exc)
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    if not cancelled and phase is None and code is not None and code != 0:
        phase = FailurePhase.EXECUTION
        message = f"Process exited with status {code}"

    durable_write_json(meta / "outcome.json", {
        "return_code": code,
        "phase": None if phase is None else phase.value,
        "message": message,
        "cancelled": cancelled,
        "published": published,
    })



def local_manifest(files, directory, identifier):
    """Derive per-file availability from immutable publication facts.

    Optional outputs that did not appear are skipped; streams must be published
    to count as available. Wildcard declarations retain their own status and
    discover concrete matches individually.
    """
    import glob
    identity = _parse_id(identifier)
    if identity is None:
        return files
    outcome = _receipt(directory, identity["run_id"], "outcome.json")
    if outcome is None:
        return files
    published = outcome.get("published", (False, False))
    for i, entry in enumerate((files.stdout, files.stderr)):
        if entry is None or entry.local_path is None:
            continue
        # When stdout/stderr are merged both declarations refer to one stream.
        success = published[i] or (
            i == 1 and files.stdout is not None and
            entry.local_path == files.stdout.local_path and published[0]
        )
        if success:
            state, message = "available", None
        elif outcome.get("phase") == FailurePhase.ARTIFACTS.value:
            state, message = "failed", outcome.get("message") or "Publication failed"
        else:
            state, message = "skipped", None
        files = files.with_file_availability(entry.local_path, state, message)
    for kind, declarations in (
        ("output", tuple(files.outputs)),
        ("cache", tuple(x.file for x in files.caches)),
    ):
        for entry in declarations:
            if entry.source_pattern is not None or entry.local_path is None:
                continue
            path = entry.local_path
            if glob.has_magic(path):
                matches = sorted(p for p in glob.glob(path, recursive=True) if os.path.isfile(p))
                for concrete in matches:
                    files = files.with_file_availability(
                        concrete, "available", kind=kind, source_pattern=path,
                    )
                files = files.with_file_availability(
                    path, "available" if matches else "skipped",
                )
            else:
                files = files.with_file_availability(
                    path, "available" if os.path.isfile(path) else "skipped",
                )
    return files

def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 3:
        raise SystemExit("usage: local_backend DIRECTORY RUN_ID JSON_SPEC")
    directory, run_id, request = args
    _worker(directory, run_id, json.loads(request))


if __name__ == "__main__":
    main()
