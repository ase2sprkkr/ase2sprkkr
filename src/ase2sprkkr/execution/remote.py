"""Remote execution through a small ASE2SPRKKR helper on another machine."""

from __future__ import annotations

import asyncio
import json
import os
import posixpath
import signal
import subprocess
import sys
import time
from pathlib import Path
from threading import Event, Lock, Thread
from uuid import uuid4

from .base import Executor, RunReference, RunState


_FINAL_STATES = {"finished", "failed", "cancelled"}
_STATE_FILE = ".ase2sprkkr-run.json"
_REQUEST_FILE = ".ase2sprkkr-request.json"


def _state_path(directory):
    return Path(directory) / _STATE_FILE


def _request_path(directory):
    return Path(directory) / _REQUEST_FILE


def _read_state(directory):
    with _state_path(directory).open() as f:
        return json.load(f)


def _write_state(run_directory, **changes):
    run_directory = Path(run_directory)
    path = _state_path(run_directory)
    try:
        with path.open() as f:
            state = json.load(f)
    except FileNotFoundError:
        state = {}
    state.update(changes)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(state))
    os.replace(temporary, path)
    return state


def submit_remote(request):
    """Start a detached direct process and return its initial state."""
    if request.get("executor", "local") != "local":
        raise ValueError("Only the local remote executor is supported for now")

    directory = Path(request["directory"]).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    request = dict(request, directory=os.fspath(directory))
    _request_path(directory).write_text(json.dumps(request))
    state = _write_state(
        directory,
        id=request["id"],
        state="queued",
        directory=os.fspath(directory),
        input_file=_absolute_file(directory, request.get("input_file")),
        output_file=_absolute_file(directory, request.get("output_file")),
        error_file=_absolute_file(directory, request.get("error_file")),
        return_code=None,
    )

    subprocess.Popen(
        [
            sys.executable,
            "-m",
            "ase2sprkkr.execution.remote_worker",
            os.fspath(directory),
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
        close_fds=True,
    )
    return state


def _absolute_file(directory, path):
    if path is None:
        return None
    path = Path(path)
    if not path.is_absolute():
        path = Path(directory) / path
    return os.fspath(path)


def status_remote(directory):
    return _read_state(Path(directory).expanduser().resolve())


def cancel_remote(directory, timeout=1.0):
    """Cancel a direct remote process group."""
    directory = Path(directory).expanduser().resolve()
    state = _read_state(directory)
    if state.get("state") in _FINAL_STATES:
        return state

    state = _write_state(directory, state="cancelled")
    pid = state.get("pid")
    if not pid:
        return state

    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return _read_state(directory)

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return _read_state(directory)
        time.sleep(0.05)

    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return _read_state(directory)


def follow_remote(directory, stdout, stderr, *, stdout_offset=0, stderr_offset=0, poll_interval=0.1):
    """Forward newly written stdout/stderr until the remote run is final."""
    directory = Path(directory).expanduser().resolve()
    offsets = {"out": stdout_offset, "err": stderr_offset}
    streams = {"out": stdout, "err": stderr}

    while True:
        state = _read_state(directory)
        had_data = False
        for kind, key in (("out", "output_file"), ("err", "error_file")):
            path = state.get(key)
            if not path:
                continue
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


def resolve_remote_mpi(mpi):
    from ..configuration import mpi_runner

    return mpi_runner(mpi)


def run_remote_worker(directory):
    """Execute one direct run inside the detached remote worker process."""
    directory = Path(directory).expanduser().resolve()
    try:
        request = json.loads(_request_path(directory).read_text())
        state = _read_state(directory)
        if state.get("state") == "cancelled":
            return

        stdin_path = request.get("stdin")
        output_path = request.get("output_file") or "stdout"
        error_path = request.get("error_file") or "stderr"
        output_path = directory / output_path
        error_path = directory / error_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        error_path.parent.mkdir(parents=True, exist_ok=True)

        stdin = open(directory / stdin_path, "rb") if stdin_path else subprocess.DEVNULL
        try:
            with output_path.open("wb") as output, error_path.open("wb") as error:
                process = subprocess.Popen(
                    request["argv"],
                    cwd=directory,
                    stdin=stdin,
                    stdout=output,
                    stderr=error,
                    start_new_session=True,
                )
                _write_state(directory, state="running", pid=process.pid)
                return_code = process.wait()
        finally:
            if stdin is not subprocess.DEVNULL:
                stdin.close()

        state = _read_state(directory)
        final_state = state.get("state")
        if final_state != "cancelled":
            final_state = "finished" if return_code == 0 else "failed"
        _write_state(directory, state=final_state, return_code=return_code)
    except BaseException as exc:
        _write_state(directory, state="failed", error=repr(exc), return_code=None)
        raise


class RemoteRunner:
    """A local proxy for one run controlled by a remote ASE2SPRKKR helper."""

    def __init__(self, executor, reader, spec, callback=None):
        self.executor = executor
        self.reader = reader
        self.spec = spec
        self.callback = callback
        self.result = reader.create_result(spec.output_file, spec.input_file)

        self._state = RunState.NEW
        self._reference = None
        self._listeners = []
        self._lock = Lock()
        self._done = Event()
        self._thread = None
        self._exception = None
        self._follow_process = None

    @property
    def reference(self):
        return self._reference

    @property
    def state(self):
        if self._reference is None or self._state in (
            RunState.NEW,
            RunState.FINISHED,
            RunState.FAILED,
            RunState.CANCELLED,
        ):
            return self._state
        status = self.executor._status(self._reference.directory)
        self._state = RunState(status.get("state", "unknown"))
        return self._state

    def add_output_listener(self, callback):
        with self._lock:
            self._listeners.append(callback)

    def _emit_output(self, data, kind):
        if not data:
            return
        text = data.decode("utf8")
        if self.reader.print_output:
            print(text, end="", file=sys.stderr if kind == "err" else sys.stdout)
        with self._lock:
            listeners = tuple(self._listeners)
        for callback in listeners:
            callback(text, kind)

    def _write_output(self, data, kind):
        target = self.spec.error_file if kind == "err" and self.spec.error_file is not None else self.spec.output_file
        if target is not None:
            target.write(data)
            target.flush()

    @staticmethod
    def _file_name(value):
        if value is None:
            return None
        if hasattr(value, "flush"):
            value.flush()
        value = getattr(value, "name", value)
        try:
            return os.fspath(value)
        except TypeError as exc:
            raise ValueError("Remote execution requires named input files") from exc

    def _stage_file(self, value, remote_directory, used):
        name = self._file_name(value)
        if name is None:
            return None
        if isinstance(name, int):
            raise ValueError("Remote execution requires named input files")
        source = os.path.abspath(name)
        target = os.path.basename(source)
        previous = used.get(target)
        if previous is not None and previous != source:
            raise ValueError("Remote input files must have distinct base names")
        used[target] = source
        self.executor.transport.put_file(source, posixpath.join(remote_directory, target))
        return target

    def _remote_argv(self, staged):
        directory = os.path.abspath(self.spec.directory)
        out = []
        for arg in self.spec.argv:
            value = os.fspath(arg)
            try:
                path = value if os.path.isabs(value) else os.path.join(directory, value)
                replacement = staged.get(os.path.abspath(path))
            except (TypeError, ValueError):
                replacement = None
            out.append(replacement or value)
        return out

    def _prepare_request(self, run_id, remote_directory):
        used = {}
        input_target = self._stage_file(self.spec.input_file, remote_directory, used)
        stdin_target = None
        staged = {source: target for target, source in used.items()}
        if self.spec.stdin is not None:
            stdin_name = self._file_name(self.spec.stdin)
            if stdin_name is not None:
                stdin_path = os.path.abspath(stdin_name)
                stdin_target = staged.get(stdin_path)
                if stdin_target is None:
                    stdin_target = self._stage_file(self.spec.stdin, remote_directory, used)
                    staged = {source: target for target, source in used.items()}

        output_name = self._file_name(self.spec.output_file)
        error_name = self._file_name(self.spec.error_file)
        output_target = os.path.basename(output_name) if output_name else "stdout"
        error_target = os.path.basename(error_name) if error_name else output_target + ".err"

        return {
            "id": run_id,
            "executor": self.executor.backend,
            "directory": remote_directory,
            "argv": self._remote_argv(staged),
            "input_file": input_target,
            "stdin": stdin_target,
            "output_file": output_target,
            "error_file": error_target,
        }

    async def _pump(self, source, target, kind):
        loop = asyncio.get_running_loop()
        try:
            while True:
                data = await loop.run_in_executor(None, source.read, 65536)
                if not data:
                    break
                self._write_output(data, kind)
                self._emit_output(data, kind)
                target.feed_data(data)
        finally:
            target.feed_eof()

    async def _monitor(self):
        stdout = asyncio.StreamReader()
        stderr = asyncio.StreamReader()
        parse_task = asyncio.create_task(self.reader.parse_streams(stdout, stderr, self.result))
        self._follow_process = self.executor._follow(self._reference.directory)
        loop = asyncio.get_running_loop()
        tasks = [
            asyncio.create_task(self._pump(self._follow_process.stdout, stdout, "out")),
            asyncio.create_task(self._pump(self._follow_process.stderr, stderr, "err")),
            loop.run_in_executor(None, self._follow_process.wait),
        ]
        try:
            _, _, follow_return_code = await asyncio.gather(*tasks)
            output, error = await parse_task
        finally:
            for task in tasks:
                if hasattr(task, "done") and not task.done():
                    task.cancel()
            if not parse_task.done():
                parse_task.cancel()
            await asyncio.gather(parse_task, return_exceptions=True)

        status = self.executor._status(self._reference.directory)
        self._state = RunState(status.get("state", "unknown"))
        if follow_return_code != 0 and self._state not in (
            RunState.FINISHED,
            RunState.FAILED,
            RunState.CANCELLED,
        ):
            raise RuntimeError("Remote output stream ended before the remote run finished")
        return_code = status.get("return_code")
        if return_code is None:
            return_code = 0 if self._state == RunState.FINISHED else 1
        return self.reader.complete_result(self.result, output, error, return_code)

    def _close_output_files(self):
        closed = set()
        for file in (self.spec.output_file, self.spec.error_file):
            if file is None or id(file) in closed:
                continue
            closed.add(id(file))
            if hasattr(file, "close"):
                file.close()

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
            if self._state not in (RunState.CANCELLED, RunState.FINISHED):
                self._state = RunState.FAILED
        finally:
            self._close_output_files()
            try:
                if self.callback:
                    self.callback(self.result)
            except BaseException as exc:
                if self._exception is None:
                    self._exception = exc
            self._done.set()

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")

            run_id = uuid4().hex
            remote_directory = posixpath.join(self.executor.root_directory, run_id)
            self.executor.transport.exec(("mkdir", "-p", "--", remote_directory))
            request = self._prepare_request(run_id, remote_directory)
            status = self.executor._submit(request)
            self._state = RunState(status["state"])
            directory = status["directory"]
            self._reference = RunReference(
                executor=f"remote:{self.executor.backend}",
                native_id=run_id,
                computer=self.executor.computer,
                directory=directory,
                input_file=status.get("input_file"),
                output_file=status.get("output_file"),
                error_file=status.get("error_file"),
            )
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def wait(self):
        if self._state == RunState.NEW:
            raise RuntimeError("Run has not been submitted")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        return self.result

    def run(self):
        if self._state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        if self._state == RunState.NEW:
            self.submit()
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.wait)

    def cancel(self):
        if self._reference is None:
            self._state = RunState.CANCELLED
            return
        status = self.executor._cancel(self._reference.directory)
        self._state = RunState(status.get("state", "cancelled"))


class RemoteExecutor(Executor):
    """Run commands through a remote ASE2SPRKKR helper and a transport."""

    name = "remote"

    def __init__(
        self,
        transport,
        *,
        backend="local",
        root_directory=".ase2sprkkr/runs",
        command=("ase2sprkkr", "remote"),
    ):
        self.transport = transport
        self.backend = backend
        self.root_directory = os.fspath(root_directory)
        self.command = tuple(command)
        self.computer = getattr(transport, "host", None)

    def _json_request(self, action, *args, payload=None):
        data = None if payload is None else json.dumps(payload).encode("utf8")
        result = self.transport.exec((*self.command, action, *args), input=data)
        return json.loads(result.stdout.decode("utf8"))

    def resolve_mpi(self, mpi=None):
        return self._json_request("resolve-mpi", payload={"mpi": mpi})["mpi"]

    def _submit(self, request):
        return self._json_request("submit", request["directory"], payload=request)

    def _status(self, directory):
        return self._json_request("status", directory)

    def _cancel(self, directory):
        return self._json_request("cancel", directory)

    def _follow(self, directory):
        return self.transport.exec_stream((*self.command, "follow", directory))

    def create_run(self, reader, spec, callback=None):
        return RemoteRunner(self, reader, spec, callback=callback)
