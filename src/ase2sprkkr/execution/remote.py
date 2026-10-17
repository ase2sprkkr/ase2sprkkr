"""Client-side remote execution staging, streaming and artifact synchronization."""

from __future__ import annotations

import asyncio
import os
import posixpath
from dataclasses import asdict, replace
from enum import Enum
from threading import Event, Lock, Thread

from .base import (
    ExecutionOutcome,
    FailurePhase,
    OutputObservers,
    RunState,
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
    """Client-side raw execution identified by an opaque backend ID.

    The remote target owns only backend facts. This runner owns staging,
    streaming and the explicit ``DOWNLOADING`` phase in which declared artifacts
    become locally available. It performs no task-specific parsing.
    """

    def __init__(
        self,
        executor,
        spec,
        callback=None,
        *,
        output_listener=None,
        start=False,
        backend_id=None,
        remote_directory=None,
        files=None,
        state=RunState.NEW,
        output_targets=None,
        cache_targets=None,
        output_patterns=None,
        cache_patterns=None,
        cache_before=None,
    ):
        self.executor = executor
        self.spec = spec
        self.callback = callback
        self._state = RunState(state)
        self._backend_id = backend_id
        self._remote_directory = remote_directory
        self._files = files
        self._output_observers = OutputObservers(output_listener)
        self._lock = Lock()
        self._done = Event()
        self._thread = None
        self._exception = None
        self._callback_exception = None
        self._return_code = None
        self._failure_phase = None
        self._failure_message = None
        self._follow_process = None
        self._output_targets = dict(output_targets or {})
        self._cache_targets = dict(cache_targets or {})
        self._output_patterns = tuple(output_patterns or ())
        self._cache_patterns = tuple(cache_patterns or ())
        self._cache_before = dict(cache_before or {})
        self.result = ExecutionOutcome(
            output_file=self._local_path(spec.output_file),
            input_file=self._local_path(spec.input_file),
            error_file=self._local_path(spec.error_file),
        )

        if start:
            self.start()

    @property
    def backend_id(self):
        return self._backend_id

    @property
    def files(self):
        return self._files

    def _set_file_availability(self, local, state, message=None, *, remote=None, kind="output"):
        if self._files is not None:
            self._files = self._files.with_file_availability(
                local, state, message, remote_path=remote, kind=kind,
            )

    def _set_pending_failed(self, message):
        """Record a transfer-wide failure on still-pending output artifacts."""
        if self._files is None:
            return
        for entry in self._files.artifacts:
            if entry.availability_state in (None, "pending"):
                self._set_file_availability(entry.local_path, "failed", message)

    @property
    def execution(self):
        if self._backend_id is None and self._remote_directory is None:
            return None
        from .endpoint import endpoint_for_executor
        from .run_model import Execution, WorkingDirectory
        return Execution(
            endpoint_for_executor(self.executor.endpoint),
            WorkingDirectory(os.path.abspath(self.spec.directory), self._remote_directory),
            backend_id=self._backend_id,
        )

    @property
    def state(self):
        # Once monitoring has started, only the client monitor may advance from
        # backend execution into DOWNLOADING. A direct backend query must never
        # expose FINISHED before artifacts are locally available.
        if self._thread is not None:
            return self._state
        if self._backend_id is None or self._state in (
            RunState.NEW,
            RunState.DOWNLOADING,
            RunState.FINISHED,
            RunState.FAILED,
            RunState.CANCELLED,
        ):
            return self._state
        status = _run_sync(
            self.executor._status(
                self._remote_directory,
                self._backend_id,
                self._files,
            )
        )
        self._state = RunState(status.get("state", "unknown"))
        if status.get("return_code") is not None:
            self._return_code = status["return_code"]
            self.result.return_code = self._return_code
        phase = status.get("failure_phase")
        self._failure_phase = None if phase is None else FailurePhase(phase)
        self._failure_message = status.get("error")
        return self._state

    @property
    def exception(self):
        return self._exception

    @property
    def return_code(self):
        return self._return_code

    @property
    def failure_phase(self):
        return self._failure_phase

    @property
    def failure_message(self):
        return self._failure_message

    def add_output_listener(self, callback):
        self._output_observers.add(callback)

    def add_output_error_listener(self, callback):
        self._output_observers.add_error_listener(callback)

    def _emit_output(self, data, kind):
        if data:
            self._output_observers.emit(data.decode("utf8"), kind)

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
        if not os.path.isabs(name):
            name = os.path.join(self.spec.directory, name)
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
                    raise FileNotFoundError(f"Local input file does not exist: {source}")
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
                        match, used, uploads=uploads, required=False
                    )
                if uploads:
                    await self.executor.transport.put_files(
                        uploads, self.executor.local_root, remote_root
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
            for path in (self.result.output_file, self.result.error_file)
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
            for entry in (() if self._files is None else self._files.artifacts):
                if entry.availability_state == "pending":
                    self._set_file_availability(entry.local_path, "skipped")
            return

        for local, remote in self._output_targets.items():
            try:
                await self._sync_one_output(local, remote, required)
            except Exception as exc:
                self._set_file_availability(local, "failed", str(exc))
                raise
            else:
                self._set_file_availability(local, "available")

        if not self._output_patterns:
            return

        remote_root = await self.executor._resolved_remote_root()
        primary = {
            path
            for path in (self.result.output_file, self.result.error_file)
            if path is not None
        }
        selected = set()
        for pattern in self._output_patterns:
            matches = await self.executor.transport.expand_files(pattern, remote_root)
            if not matches:
                if required:
                    self._set_file_availability(
                        os.fspath(await self.executor.unmap_path(pattern)), "failed",
                        f"Declared remote output pattern did not match any files: {pattern}",
                        remote=pattern,
                    )
                    raise FileNotFoundError(
                        f"Declared remote output pattern did not match any files: {pattern}"
                    )
                continue
            for remote in matches:
                if remote in selected:
                    continue
                selected.add(remote)
                local = os.fspath(await self.executor.unmap_path(remote))
                if local in primary:
                    continue
                try:
                    await self._sync_one_output(local, remote, required)
                except Exception as exc:
                    self._set_file_availability(
                        local, "failed", str(exc), remote=remote,
                    )
                    raise
                else:
                    self._set_file_availability(
                        local, "available", remote=remote,
                    )

    async def _sync_one_output(self, local, remote, required):
        policy = self.executor.output_files
        remote_digest = await self.executor.transport.file_digest(remote)
        if remote_digest is None:
            if required:
                raise FileNotFoundError(
                    f"Declared remote output file was not created: {remote}"
                )
            return
        if policy == FileSyncPolicy.ALWAYS:
            await self.executor.transport.get_file(remote, local)
            return
        local_digest = local_file_digest(local)
        if policy == FileSyncPolicy.MISSING:
            if local_digest is None:
                await self.executor.transport.get_file(remote, local)
        elif local_digest != remote_digest:
            await self.executor.transport.get_file(remote, local)

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
            for remote in await self.executor.transport.expand_files(pattern, remote_root):
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
                downloads, remote_root, self.executor.local_root
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

        input_path = self.result.input_file
        input_target = staged.get(input_path)
        if input_target is None:
            input_target = await self._stage_file(self.spec.input_file, used)
            staged = {source: target for target, source in used.items()}

        stdin_target = None
        if self.spec.stdin is not None:
            stdin_path = self._local_path(self.spec.stdin)
            stdin_target = staged.get(stdin_path)
            if stdin_target is None:
                stdin_target = await self._stage_file(self.spec.stdin, used)
                staged = {source: target for target, source in used.items()}

        await self._prepare_cache_files(used)
        staged = {source: target for target, source in used.items()}
        await self._prepare_output_files(set(used))

        output_target = (
            await self.executor.map_path(self.result.output_file)
            if self.result.output_file
            else posixpath.join(remote_directory, "stdout")
        )
        error_target = (
            await self.executor.map_path(self.result.error_file)
            if self.result.error_file
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
            "mpi_argv": (
                None
                if self.spec.mpi_argv is None
                else self._remote_argv(self.spec.mpi_argv, staged)
            ),
            "launch": asdict(self.spec.launch),
            "input_file": input_target,
            "stdin": stdin_target,
            "output_file": output_target,
            "error_file": error_target,
            "client_sync": client_sync,
        }

    async def _pump(self, source, kind):
        while True:
            data = await source.read(65536)
            if not data:
                return
            self._write_output(data, kind)
            self._emit_output(data, kind)

    async def _monitor(self):
        stdout_task = stderr_task = None
        outputs_committed = False
        try:
            self._follow_process = await self.executor._follow(
                self._remote_directory,
                self._backend_id,
                self._files,
            )
            stdout_task = asyncio.create_task(
                self._pump(self._follow_process.stdout, "out")
            )
            stderr_task = asyncio.create_task(
                self._pump(self._follow_process.stderr, "err")
            )
            follow_return_code, _, _ = await asyncio.gather(
                self._follow_process.wait(), stdout_task, stderr_task
            )

            status = await self.executor._status(
                self._remote_directory,
                self._backend_id,
                self._files,
            )
            backend_state = RunState(status.get("state", "unknown"))
            self._return_code = status.get("return_code")
            self.result.return_code = self._return_code
            phase = status.get("failure_phase")
            self._failure_phase = None if phase is None else FailurePhase(phase)
            self._failure_message = status.get("error")

            if follow_return_code != 0 and backend_state not in (
                RunState.FINISHED,
                RunState.FAILED,
                RunState.CANCELLED,
            ):
                raise RuntimeError(
                    "Remote output stream ended before the remote run finished"
                )

            if backend_state == RunState.FAILED and self._failure_phase in (
                FailurePhase.SUBMISSION,
                FailurePhase.LAUNCH,
            ):
                self._state = RunState.FAILED
                return self.result

            if backend_state not in (
                RunState.FINISHED,
                RunState.FAILED,
                RunState.CANCELLED,
            ):
                # The stream helper may stop because a backend identity was lost
                # (for example a direct-process hard crash before its durable
                # outcome barrier). Without an authoritative terminal backend
                # fact, streamed bytes are provisional and must not replace the
                # previously published local artifacts.
                self._state = RunState.UNKNOWN
                self._failure_phase = None
                if self._failure_message is None:
                    self._failure_message = (
                        "Remote output stream ended without a terminal backend "
                        f"state ({backend_state.value})"
                    )
                return self.result

            self._state = RunState.DOWNLOADING
            try:
                await self._sync_output_files(required=backend_state == RunState.FINISHED)
                await self._sync_cache_files()
                self._finish_output_files(commit=True)
            except Exception as exc:
                self._set_pending_failed(str(exc))
                if backend_state == RunState.FINISHED:
                    # A successful calculation with unavailable artifacts is a
                    # failed run, even though its backend return code is zero.
                    self._state = RunState.FAILED
                    self._failure_phase = FailurePhase.DOWNLOAD
                    self._failure_message = str(exc)
                    raise
                # Diagnostic download failure cannot replace the original
                # failed/cancelled backend outcome or its failure reason.
                self._state = backend_state
                if backend_state == RunState.FAILED and self._failure_phase is None:
                    self._failure_phase = FailurePhase.EXECUTION
                return self.result

            outputs_committed = True
            if self._files is not None:
                for entry in self._files.artifacts:
                    if entry.availability_state == "pending":
                        self._set_file_availability(entry.local_path, "available")
            self._state = backend_state
            if self._state == RunState.FAILED and self._failure_phase is None:
                self._failure_phase = FailurePhase.EXECUTION
            return self.result
        except BaseException as exc:
            if self._failure_phase is None:
                self._failure_phase = FailurePhase.DOWNLOAD
            if self._failure_message is None:
                self._failure_message = str(exc)
            if self._state != RunState.CANCELLED:
                self._state = RunState.FAILED
            raise
        finally:
            tasks = tuple(
                task for task in (stdout_task, stderr_task) if task is not None
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

    def _thread_main(self):
        try:
            asyncio.run(self._monitor())
        except BaseException as exc:
            self._exception = exc
        finally:
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
            if self._backend_id is None:
                raise RuntimeError("Run has no backend reference")
            self._thread = Thread(target=self._thread_main, daemon=True)
            self._thread.start()

    def start(self):
        if self._backend_id is None:
            return self.submit()
        return self._start_monitor()

    def submit(self):
        raise NotImplementedError(
            "Remote submission is provided by the endpoint-specific runner"
        )

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
        return await asyncio.to_thread(self.run)

    def cancel(self):
        if self._backend_id is None:
            raise RuntimeError("Run has not been submitted")
        status = _run_sync(
            self.executor._cancel(
                self._remote_directory,
                self._backend_id,
            )
        )
        self._state = RunState(status.get("state", "cancelled"))
        if status.get("return_code") is not None:
            self._return_code = status["return_code"]
            self.result.return_code = self._return_code
        phase = status.get("failure_phase")
        self._failure_phase = None if phase is None else FailurePhase(phase)
        self._failure_message = status.get("error")
