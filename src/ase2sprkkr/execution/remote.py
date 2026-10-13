"""Remote execution through a small ASE2SPRKKR helper on another machine."""

from __future__ import annotations

import asyncio
import os
import posixpath
import sys
from dataclasses import asdict
from enum import Enum
from threading import Event, Lock, Thread
from uuid import uuid4

from .base import (
    OutputObservers,
    RunReference,
    RunState,
    SubmissionNotStarted,
)
from .detached import (
    _FINAL_STATES,
    _METADATA_DIR,
    _REQUEST_FILE,
    _RUNS_DIR,
    _STATE_FILE,
    _STATE_LOCK_FILE,
    _DETACHED_EXECUTORS as _REMOTE_EXECUTORS,
    _read_state,
    _read_state_unlocked,
    _request_path,
    _request_spec,
    _run_directory,
    _state_lock,
    _state_lock_path,
    _state_path,
    _timestamp,
    _write_state,
    _write_state_unlocked,
    cancel_detached as cancel_remote,
    clean_detached_runs as clean_remote_runs,
    follow_detached as follow_remote,
    list_detached_runs as list_remote_runs,
    resolve_detached_mpi as resolve_remote_mpi,
    run_detached_worker as run_remote_worker,
    status_detached as status_remote,
    submit_detached as submit_remote,
)
from .transport import expand_local_files, has_file_magic, local_file_digest


class FileSyncPolicy(str, Enum):
    """Policy controlling transfers between local and remote workspaces."""

    CHANGED = "changed"
    ALWAYS = "always"
    MISSING = "missing"
    NEVER = "never"


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

    def __init__(
        self,
        executor,
        reader,
        spec,
        callback=None,
        *,
        output_listener=None,
        start=False,
        reference=None,
        state=RunState.NEW,
        output_targets=None,
        cache_targets=None,
        output_patterns=None,
        cache_patterns=None,
        cache_before=None,
    ):
        self.executor = executor
        self.reader = reader
        self.spec = spec
        self.callback = callback
        self.result = reader.create_result(spec.output_file, spec.input_file)

        self._state = RunState(state)
        self._reference = reference
        self._output_observers = OutputObservers(output_listener)
        self._lock = Lock()
        self._done = Event()
        self._thread = None
        self._exception = None
        self._callback_exception = None
        self._return_code = None
        self._file_retrieval_state = "pending"
        self._file_retrieval_message = None
        self._follow_process = None
        self._output_targets = dict(output_targets or {})
        self._cache_targets = dict(cache_targets or {})
        self._output_patterns = tuple(output_patterns or ())
        self._cache_patterns = tuple(cache_patterns or ())
        self._cache_before = dict(cache_before or {})

        if start:
            self.start()

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
        status = _run_sync(
            self.executor._status(
                self._reference.directory, self._reference.native_id
            )
        )
        self._state = RunState(status.get("state", "unknown"))
        if status.get("return_code") is not None:
            self._return_code = status["return_code"]
        return self._state

    @property
    def return_code(self):
        return self._return_code

    @property
    def file_retrieval_state(self):
        return self._file_retrieval_state

    @property
    def file_retrieval_message(self):
        return self._file_retrieval_message

    def add_output_listener(self, callback):
        self._output_observers.add(callback)

    def add_output_error_listener(self, callback):
        self._output_observers.add_error_listener(callback)

    def _emit_output(self, data, kind):
        if not data:
            return
        text = data.decode("utf8")
        if self.reader.print_output:
            print(text, end="", file=sys.stderr if kind == "err" else sys.stdout)
        self._output_observers.emit(text, kind)

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
            raise ValueError("Remote execution requires named files") from exc

    def _local_path(self, value):
        name = self._file_name(value)
        if name is None:
            return None
        if isinstance(name, int):
            raise ValueError("Remote execution requires named files")
        return os.path.abspath(name)

    async def _stage_file(self, value, used, *, uploads=None, required=True):
        source = self._local_path(value)
        if source is None:
            return None
        target = await self.executor.map_path(source)
        previous = used.get(target)
        if previous is not None and previous != source:
            raise ValueError(f"Remote input target collision: {target}")
        if previous is not None:
            return target

        policy = self.executor.input_files
        if policy == FileSyncPolicy.NEVER:
            used[target] = source
            return target

        local_digest = local_file_digest(source)
        remote_digest = await self.executor.transport.file_digest(target)
        transfer = False

        if policy == FileSyncPolicy.ALWAYS:
            if local_digest is None:
                if required:
                    raise FileNotFoundError(
                        f"Local input file does not exist: {source}"
                    )
            else:
                transfer = True
        elif policy == FileSyncPolicy.MISSING:
            if remote_digest is None:
                if local_digest is None:
                    if required:
                        raise FileNotFoundError(
                            f"Input file exists neither locally nor remotely: {source}"
                        )
                else:
                    transfer = True
        elif local_digest is None:
            if remote_digest is None and required:
                raise FileNotFoundError(
                    f"Input file exists neither locally nor remotely: {source}"
                )
        elif local_digest != remote_digest:
            transfer = True

        if transfer:
            if uploads is None:
                await self.executor.transport.put_file(source, target)
            else:
                uploads.append(source)

        used[target] = source
        return target

    async def _stage_input_pattern(self, value, used):
        source = self._local_path(value)
        matches = expand_local_files(source, self.executor.local_root)
        if not matches:
            raise FileNotFoundError(
                f"Local input pattern did not match any files: {source}"
            )

        uploads = []
        for match in matches:
            await self._stage_file(match, used, uploads=uploads)
        if uploads:
            await self.executor.transport.put_files(
                uploads,
                self.executor.local_root,
                await self.executor._resolved_remote_root(),
            )

    async def _prepare_cache_files(self, used):
        targets = {}
        patterns = []
        pattern_before = {}
        input_policy = self.executor.input_files
        output_policy = self.executor.output_files
        remote_root = await self.executor._resolved_remote_root()

        for value in self.spec.cache_files:
            local = self._local_path(value)
            if local is None:
                continue

            if has_file_magic(local):
                remote_pattern = await self.executor.map_path(local)
                if remote_pattern in patterns:
                    continue
                patterns.append(remote_pattern)

                uploads = []
                for match in expand_local_files(local, self.executor.local_root):
                    await self._stage_file(
                        match,
                        used,
                        uploads=uploads,
                        required=False,
                    )
                if uploads:
                    await self.executor.transport.put_files(
                        uploads,
                        self.executor.local_root,
                        remote_root,
                    )

                if output_policy != FileSyncPolicy.NEVER:
                    for remote in await self.executor.transport.expand_files(
                        remote_pattern, remote_root
                    ):
                        pattern_before[remote] = (
                            await self.executor.transport.file_signature(remote)
                        )
                continue

            if local in targets:
                continue
            remote = await self.executor.map_path(local)
            previous = used.get(remote)
            if previous is not None and previous != local:
                raise ValueError(f"Remote cache target collision: {remote}")

            if previous is None:
                if input_policy != FileSyncPolicy.NEVER:
                    local_digest = local_file_digest(local)
                    remote_digest = await self.executor.transport.file_digest(remote)
                    if input_policy == FileSyncPolicy.ALWAYS:
                        if local_digest is not None:
                            await self.executor.transport.put_file(local, remote)
                    elif input_policy == FileSyncPolicy.MISSING:
                        if remote_digest is None and local_digest is not None:
                            await self.executor.transport.put_file(local, remote)
                    elif local_digest is not None and local_digest != remote_digest:
                        await self.executor.transport.put_file(local, remote)
                used[remote] = local

            before = None
            if output_policy != FileSyncPolicy.NEVER:
                before = await self.executor.transport.file_signature(remote)
            targets[local] = (remote, before)

        self._cache_targets = targets
        self._cache_patterns = tuple(patterns)
        self._cache_before = pattern_before

    async def _prepare_output_files(self, input_targets):
        primary = {
            path
            for path in (
                self._local_path(self.spec.output_file),
                self._local_path(self.spec.error_file),
            )
            if path is not None
        }
        targets = {}
        patterns = []
        remote_root = await self.executor._resolved_remote_root()
        removable = set()

        for value in self.spec.output_files:
            local = self._local_path(value)
            if local is None or local in primary:
                continue

            remote = await self.executor.map_path(local)
            if has_file_magic(local):
                if remote in patterns:
                    continue
                patterns.append(remote)
                for match in await self.executor.transport.expand_files(
                    remote, remote_root
                ):
                    if match not in input_targets:
                        removable.add(match)
                continue

            if local in targets:
                continue
            targets[local] = remote
            if remote not in input_targets:
                await self.executor.transport.run(("rm", "-f", "--", remote))

        if removable:
            await self.executor.transport.run(
                ("rm", "-f", "--", *sorted(removable))
            )
        self._output_targets = targets
        self._output_patterns = tuple(patterns)

    async def _sync_output_files(self, required):
        policy = self.executor.output_files
        if policy == FileSyncPolicy.NEVER:
            return

        for local, remote in self._output_targets.items():
            remote_digest = await self.executor.transport.file_digest(remote)
            if remote_digest is None:
                if required:
                    raise FileNotFoundError(
                        f"Declared remote output file was not created: {remote}"
                    )
                continue

            if policy == FileSyncPolicy.ALWAYS:
                await self.executor.transport.get_file(remote, local)
                continue
            local_digest = local_file_digest(local)
            if policy == FileSyncPolicy.MISSING:
                if local_digest is None:
                    await self.executor.transport.get_file(remote, local)
            elif local_digest != remote_digest:
                await self.executor.transport.get_file(remote, local)

        if not self._output_patterns:
            return

        remote_root = await self.executor._resolved_remote_root()
        primary = {
            path
            for path in (
                self._local_path(self.spec.output_file),
                self._local_path(self.spec.error_file),
            )
            if path is not None
        }
        downloads = []
        selected = set()
        for pattern in self._output_patterns:
            matches = await self.executor.transport.expand_files(pattern, remote_root)
            if not matches:
                if required:
                    raise FileNotFoundError(
                        f"Declared remote output pattern did not match any files: {pattern}"
                    )
                continue

            for remote in matches:
                if remote in selected:
                    continue
                local = os.fspath(await self.executor.unmap_path(remote))
                if local in primary:
                    selected.add(remote)
                    continue

                download = policy == FileSyncPolicy.ALWAYS
                if not download:
                    local_digest = local_file_digest(local)
                    if policy == FileSyncPolicy.MISSING:
                        download = local_digest is None
                    else:
                        remote_digest = await self.executor.transport.file_digest(remote)
                        download = local_digest != remote_digest
                if download:
                    downloads.append(remote)
                selected.add(remote)

        if downloads:
            await self.executor.transport.get_files(
                downloads,
                remote_root,
                self.executor.local_root,
            )

    async def _sync_cache_files(self):
        policy = self.executor.output_files
        if policy == FileSyncPolicy.NEVER:
            return

        for local, (remote, before) in self._cache_targets.items():
            after = await self.executor.transport.file_signature(remote)
            if after is None or after == before:
                continue
            if policy == FileSyncPolicy.MISSING and local_file_digest(local) is not None:
                continue
            await self.executor.transport.get_file(remote, local)

        if not self._cache_patterns:
            return

        remote_root = await self.executor._resolved_remote_root()
        downloads = []
        selected = set()
        exact = {remote for remote, _ in self._cache_targets.values()}
        for pattern in self._cache_patterns:
            for remote in await self.executor.transport.expand_files(
                pattern, remote_root
            ):
                if remote in selected or remote in exact:
                    continue
                selected.add(remote)
                after = await self.executor.transport.file_signature(remote)
                if after is None or after == self._cache_before.get(remote):
                    continue
                local = os.fspath(await self.executor.unmap_path(remote))
                if policy == FileSyncPolicy.MISSING and local_file_digest(local) is not None:
                    continue
                downloads.append(remote)

        if downloads:
            await self.executor.transport.get_files(
                downloads,
                remote_root,
                self.executor.local_root,
            )

    def _remote_argv(self, argv, staged):
        directory = os.path.abspath(self.spec.directory)
        out = []
        for arg in argv:
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
            source = self._local_path(value)
            if source is not None and has_file_magic(source):
                await self._stage_input_pattern(value, used)
            else:
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

        await self._prepare_cache_files(used)
        staged = {source: target for target, source in used.items()}
        await self._prepare_output_files(set(used))
        output_name = self._file_name(self.spec.output_file)
        error_name = self._file_name(self.spec.error_file)
        output_target = (
            await self.executor.map_path(output_name)
            if output_name
            else posixpath.join(remote_directory, "stdout")
        )
        error_target = (
            await self.executor.map_path(error_name)
            if error_name
            else output_target + ".err"
        )
        client_sync = {
            "output_file": output_target if self.spec.output_file is not None else None,
            "error_file": error_target if self.spec.error_file is not None else None,
            "output_files": list(self._output_targets.values()),
            "cache_files": [
                {"path": remote, "before": before}
                for remote, before in self._cache_targets.values()
            ],
            "output_patterns": list(self._output_patterns),
            "cache_patterns": list(self._cache_patterns),
            "cache_before": [
                {"path": remote, "signature": before}
                for remote, before in self._cache_before.items()
            ],
        }

        return {
            "id": run_id,
            "executor": self.executor.backend,
            "executor_options": self.executor.backend_options,
            "directory": remote_directory,
            "argv": self._remote_argv(self.spec.argv, staged),
            "mpi_argv": None
            if self.spec.mpi_argv is None
            else self._remote_argv(self.spec.mpi_argv, staged),
            "launch": asdict(self.spec.launch),
            "input_file": input_target,
            "stdin": stdin_target,
            "output_file": output_target,
            "error_file": error_target,
            "client_sync": client_sync,
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
        stdout_task = stderr_task = None
        outputs_committed = False
        try:
            try:
                self._follow_process = await self.executor._follow(
                    self._reference.directory, self._reference.native_id
                )
                stdout_task = asyncio.create_task(
                    self._pump(self._follow_process.stdout, stdout, "out")
                )
                stderr_task = asyncio.create_task(
                    self._pump(self._follow_process.stderr, stderr, "err")
                )
                follow_return_code, _, _ = await asyncio.gather(
                    self._follow_process.wait(), stdout_task, stderr_task
                )
            except Exception as exc:
                self._file_retrieval_state = "failed"
                self._file_retrieval_message = str(exc)
                try:
                    status = await self.executor._status(
                        self._reference.directory, self._reference.native_id
                    )
                    self._state = RunState(status.get("state", "unknown"))
                    if status.get("return_code") is not None:
                        self._return_code = status["return_code"]
                except Exception:
                    pass
                raise

            status = await self.executor._status(
                self._reference.directory, self._reference.native_id
            )
            self._state = RunState(status.get("state", "unknown"))
            if follow_return_code != 0 and self._state not in (
                RunState.FINISHED,
                RunState.FAILED,
                RunState.CANCELLED,
            ):
                exc = RuntimeError(
                    "Remote output stream ended before the remote run finished"
                )
                self._file_retrieval_state = "failed"
                self._file_retrieval_message = str(exc)
                raise exc
            self._return_code = status.get("return_code")
            if self._return_code is None:
                self._return_code = 0 if self._state == RunState.FINISHED else 1

            try:
                await self._sync_output_files(required=self._state == RunState.FINISHED)
                await self._sync_cache_files()
            except Exception as exc:
                self._file_retrieval_state = "failed"
                self._file_retrieval_message = str(exc)
                raise
            else:
                self._file_retrieval_state = "available"
                self._file_retrieval_message = None

            self._finish_output_files(commit=True)
            outputs_committed = True

            output, error = await parse_task
            return self.reader.complete_result(
                self.result, output, error, self._return_code
            )
        finally:
            tasks = tuple(
                task
                for task in (stdout_task, stderr_task, parse_task)
                if task is not None
            )
            for task in tasks:
                if not task.done():
                    task.cancel()
            if tasks:
                await asyncio.gather(*tasks, return_exceptions=True)
            if not outputs_committed:
                self._finish_output_files(commit=False)

    def _finish_output_files(self, *, commit):
        closed = set()
        for file in (self.spec.output_file, self.spec.error_file):
            if file is None or id(file) in closed:
                continue
            closed.add(id(file))
            method = getattr(file, "commit" if commit else "discard", None)
            if method is not None:
                method()
            elif hasattr(file, "close"):
                file.close()

    def _close_output_files(self):
        self._finish_output_files(commit=False)

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
        finally:
            self._close_output_files()
            try:
                if self.callback:
                    self.callback(self.result)
            except BaseException as exc:
                self._callback_exception = exc
            self._done.set()

    def _start_monitor(self):
        with self._lock:
            if self._thread is not None:
                return
            if self._reference is None:
                raise RuntimeError("Run has no backend reference")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def start(self):
        """Start submission or monitoring, depending on runner state."""
        if self._reference is None:
            return self.submit()
        return self._start_monitor()

    async def _submit_async(self):
        run_id = uuid4().hex
        try:
            remote_directory = await self.executor.map_directory(self.spec.directory)
            await self.executor.transport.run(
                ("mkdir", "-p", "--", remote_directory)
            )
            request = await self._prepare_request(run_id, remote_directory)
        except Exception as exc:
            raise SubmissionNotStarted(str(exc)) from exc
        status = await self.executor._submit(request)
        return run_id, status

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")

            run_id, status = _run_sync(self._submit_async())
            self._state = RunState(status["state"])
            self._return_code = status.get("return_code")
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
        self._start_monitor()

    def wait(self):
        if self._thread is None:
            if self._state == RunState.NEW:
                raise RuntimeError("Run has not been submitted")
            raise RuntimeError("Attached run has not been started")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        if self._callback_exception is not None:
            raise self._callback_exception
        return self.result

    def run(self):
        if self._thread is None:
            self.start()
        return self.wait()

    async def run_async(self):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.run)

    def cancel(self):
        if self._reference is None:
            raise RuntimeError("Run has not been submitted")
        status = _run_sync(
            self.executor._cancel(
                self._reference.directory, self._reference.native_id
            )
        )
        self._state = RunState(status.get("state", "cancelled"))
        if status.get("return_code") is not None:
            self._return_code = status["return_code"]
