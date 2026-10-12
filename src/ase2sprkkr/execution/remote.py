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
    try:
        previous = _read_state(directory)
    except FileNotFoundError:
        previous = None
    if previous and previous.get("state") not in _FINAL_STATES:
        raise RuntimeError(
            f"Remote workspace already contains an active run: {directory}"
        )

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

    _request_path(directory).write_text(json.dumps(request))
    state = _write_state(
        directory,
        id=request["id"],
        state="queued",
        directory=os.fspath(directory),
        input_file=_absolute_file(directory, request.get("input_file")),
        output_file=output_file,
        error_file=error_file,
        return_code=None,
        error=None,
        pid=None,
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


def follow_remote(
    directory,
    stdout,
    stderr,
    *,
    stdout_offset=0,
    stderr_offset=0,
    poll_interval=0.1,
):
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


def _run_sync(coro):
    """Run one transport coroutine from the synchronous executor facade."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)

    result = []
    error = []

    def worker():
        try:
            result.append(asyncio.run(coro))
        except BaseException as exc:
            error.append(exc)

    thread = Thread(target=worker)
    thread.start()
    thread.join()
    if error:
        raise error[0]
    return result[0]


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
        status = _run_sync(self.executor._status(self._reference.directory))
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
        target = (
            self.spec.error_file
            if kind == "err" and self.spec.error_file is not None
            else self.spec.output_file
        )
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

    async def _stage_file(self, value, used):
        name = self._file_name(value)
        if name is None:
            return None
        if isinstance(name, int):
            raise ValueError("Remote execution requires named input files")
        source = os.path.abspath(name)
        target = await self.executor.map_path(source)
        previous = used.get(target)
        if previous is not None and previous != source:
            raise ValueError(f"Remote input target collision: {target}")
        if previous is None:
            used[target] = source
            await self.executor.transport.put_file(source, target)
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

    async def _prepare_request(self, run_id, remote_directory):
        used = {}
        for value in self.spec.input_files:
            await self._stage_file(value, used)
        staged = {source: target for target, source in used.items()}

        input_name = self._file_name(self.spec.input_file)
        input_path = os.path.abspath(input_name) if input_name is not None else None
        input_target = staged.get(input_path)
        if input_target is None:
            input_target = await self._stage_file(self.spec.input_file, used)
            staged = {source: target for target, source in used.items()}

        stdin_target = None
        if self.spec.stdin is not None:
            stdin_name = self._file_name(self.spec.stdin)
            if stdin_name is not None:
                stdin_path = os.path.abspath(stdin_name)
                stdin_target = staged.get(stdin_path)
                if stdin_target is None:
                    stdin_target = await self._stage_file(self.spec.stdin, used)
                    staged = {source: target for target, source in used.items()}

        output_name = self._file_name(self.spec.output_file)
        error_name = self._file_name(self.spec.error_file)
        output_target = (
            await self.executor.map_path(output_name) if output_name else posixpath.join(remote_directory, "stdout")
        )
        error_target = (
            await self.executor.map_path(error_name)
            if error_name
            else output_target + ".err"
        )

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
        try:
            while True:
                data = await source.read(65536)
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
        parse_task = asyncio.create_task(
            self.reader.parse_streams(stdout, stderr, self.result)
        )
        self._follow_process = await self.executor._follow(
            self._reference.directory
        )
        stdout_task = asyncio.create_task(
            self._pump(self._follow_process.stdout, stdout, "out")
        )
        stderr_task = asyncio.create_task(
            self._pump(self._follow_process.stderr, stderr, "err")
        )
        try:
            follow_return_code, _, _ = await asyncio.gather(
                self._follow_process.wait(), stdout_task, stderr_task
            )
            output, error = await parse_task
        finally:
            for task in (stdout_task, stderr_task, parse_task):
                if not task.done():
                    task.cancel()
            await asyncio.gather(
                stdout_task, stderr_task, parse_task, return_exceptions=True
            )

        status = await self.executor._status(self._reference.directory)
        self._state = RunState(status.get("state", "unknown"))
        if follow_return_code != 0 and self._state not in (
            RunState.FINISHED,
            RunState.FAILED,
            RunState.CANCELLED,
        ):
            raise RuntimeError(
                "Remote output stream ended before the remote run finished"
            )
        return_code = status.get("return_code")
        if return_code is None:
            return_code = 0 if self._state == RunState.FINISHED else 1
        return self.reader.complete_result(
            self.result, output, error, return_code
        )

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

    async def _submit_async(self):
        run_id = uuid4().hex
        remote_directory = await self.executor.map_directory(self.spec.directory)
        await self.executor.transport.run(
            ("mkdir", "-p", "--", remote_directory)
        )
        request = await self._prepare_request(run_id, remote_directory)
        status = await self.executor._submit(request)
        return run_id, status

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")

            run_id, status = _run_sync(self._submit_async())
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
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.run)

    def cancel(self):
        if self._reference is None:
            self._state = RunState.CANCELLED
            return
        status = _run_sync(self.executor._cancel(self._reference.directory))
        self._state = RunState(status.get("state", "cancelled"))


class RemoteExecutor(Executor):
    """Run commands through a remote ASE2SPRKKR helper and a transport.

    Calculation directories are persistent workspaces.  By default a local
    path below the local home directory is mapped to the same relative path
    below the remote home directory.  ``local_root`` and ``remote_root`` can
    override either side of that mapping.
    """

    name = "remote"

    def __init__(
        self,
        transport,
        *,
        backend="local",
        local_root=None,
        remote_root=None,
        root_directory=None,
        command=("ase2sprkkr", "remote"),
    ):
        if root_directory is not None:
            if remote_root is not None:
                raise ValueError("Use either remote_root or root_directory, not both")
            remote_root = root_directory
        self.transport = transport
        self.backend = backend
        self.local_root = Path(local_root or Path.home()).expanduser().resolve()
        self.remote_root = None if remote_root is None else os.fspath(remote_root)
        self.command = tuple(command)
        self.computer = getattr(transport, "host", None)
        self._remote_home_cache = None

    async def _remote_home(self):
        if self._remote_home_cache is None:
            result = await self.transport.run(
                ("sh", "-c", 'printf "%s" "$HOME"')
            )
            home = result.stdout.decode("utf8").strip()
            if not home:
                raise RuntimeError("Could not determine remote home directory")
            self._remote_home_cache = posixpath.normpath(home)
        return self._remote_home_cache

    async def _resolved_remote_root(self):
        home = await self._remote_home()
        if self.remote_root is None:
            return home
        root = self.remote_root
        if root == "~":
            return home
        if root.startswith("~/"):
            return posixpath.join(home, root[2:])
        if not posixpath.isabs(root):
            return posixpath.join(home, root)
        return posixpath.normpath(root)

    async def map_path(self, path):
        """Map one local path to the corresponding path on the target filesystem."""
        path = Path(path).expanduser().resolve()
        try:
            relative = path.relative_to(self.local_root)
        except ValueError as exc:
            raise ValueError(
                f"Local path {path} is outside mapped root {self.local_root}; "
                "specify local_root/remote_root explicitly"
            ) from exc
        remote_root = await self._resolved_remote_root()
        if not relative.parts:
            return remote_root
        return posixpath.join(remote_root, *relative.parts)

    async def map_directory(self, directory):
        """Map a local calculation directory to its persistent remote workspace."""
        return await self.map_path(directory)

    async def _json_request(self, action, *args, payload=None):
        data = None if payload is None else json.dumps(payload).encode("utf8")
        result = await self.transport.run(
            (*self.command, action, *args), input=data
        )
        return json.loads(result.stdout.decode("utf8"))

    def resolve_mpi(self, mpi=None):
        return _run_sync(
            self._json_request("resolve-mpi", payload={"mpi": mpi})
        )["mpi"]

    async def _submit(self, request):
        return await self._json_request(
            "submit", request["directory"], payload=request
        )

    async def _status(self, directory):
        return await self._json_request("status", directory)

    async def _cancel(self, directory):
        return await self._json_request("cancel", directory)

    async def _follow(self, directory):
        return await self.transport.start(
            (*self.command, "follow", directory)
        )

    def create_run(self, reader, spec, callback=None):
        return RemoteRunner(self, reader, spec, callback=callback)
