"""Reusable remote endpoint with per-run staging configuration."""

from __future__ import annotations

import asyncio
import glob
import json
import os
import posixpath
import warnings
from pathlib import Path, PurePosixPath
from uuid import uuid4

from .base import (
    Executor,
    LaunchSpec,
    MpiNotFoundWarning,
    RunReference,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .remote import (
    FileSyncPolicy,
    RemoteRunner,
    _METADATA_DIR,
    _REQUEST_FILE,
    _REMOTE_EXECUTORS,
    _RUNS_DIR,
    _STATE_FILE,
    _run_sync,
)
from .transport import AtomicLocalFile


class _RemoteRun(RemoteRunner):
    """Remote runner that can reuse an identity chosen before submission."""

    def __init__(self, *args, run_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._run_id = run_id
        self._streamed_sizes = {"out": 0, "err": 0}

    def _write_output(self, data, kind):
        super()._write_output(data, kind)
        self._streamed_sizes[kind] += len(data)

    async def _monitor(self):
        # Resolve the local path while the stream target is still open. The base
        # monitor commits/closes it before returning.
        local_output = self._local_path(self.spec.output_file)
        result = await super()._monitor()
        reference = self.reference
        if reference is None:
            return result

        signatures = {}
        for kind, remote in (
            ("out", reference.output_file),
            ("err", reference.error_file),
        ):
            if remote is None:
                continue
            signature = await self.executor.transport.file_signature(remote)
            if signature is None:
                continue
            signatures[kind] = signature
            if self._streamed_sizes[kind] != signature[0]:
                stream = "stdout" if kind == "out" else "stderr"
                raise RuntimeError(
                    f"Streamed {stream} size does not match the authoritative "
                    f"remote {stream}: {remote}"
                )

        # The public output file may intentionally contain stderr as well when
        # no separate error_file was requested, so its size is not a valid
        # stdout integrity check. Preserve the authoritative stdout mtime only.
        if local_output is not None and "out" in signatures:
            stat = os.stat(local_output)
            os.utime(local_output, ns=(stat.st_atime_ns, signatures["out"][1]))
        return result

    async def _submit_async(self):
        run_id = self._run_id or uuid4().hex
        try:
            remote_directory = await self.executor.map_directory(self.spec.directory)
            await self.executor.transport.run(("mkdir", "-p", "--", remote_directory))
            request = await self._prepare_request(run_id, remote_directory)
        except Exception as exc:
            raise SubmissionNotStarted(str(exc)) from exc
        status = await self.executor._submit(request)
        return run_id, status


class _RemoteRunContext:
    """Per-run filesystem/staging policy for one reusable remote endpoint."""

    def __init__(self, endpoint, launch):
        self.endpoint = endpoint
        self.transport = endpoint.transport
        self.backend = endpoint.backend
        self.backend_options = endpoint.backend_options
        self.command = endpoint.command
        self.computer = endpoint.computer
        self.local_root = Path(launch.local_root or Path.home()).expanduser().resolve()
        self.remote_root = None if launch.remote_root is None else os.fspath(launch.remote_root)
        self.input_files = FileSyncPolicy(launch.input_sync)
        self.output_files = FileSyncPolicy(launch.output_sync)

    async def _remote_home(self):
        return await self.endpoint._remote_home()

    async def _resolved_remote_root(self):
        home = await self._remote_home()
        if self.remote_root is None or self.remote_root == "~":
            return home
        if self.remote_root.startswith("~/"):
            return posixpath.join(home, self.remote_root[2:])
        if not posixpath.isabs(self.remote_root):
            return posixpath.join(home, self.remote_root)
        return posixpath.normpath(self.remote_root)

    async def map_path(self, path):
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

    async def unmap_path(self, path):
        remote_root = await self._resolved_remote_root()
        path = posixpath.normpath(os.fspath(path))
        if not posixpath.isabs(path):
            raise ValueError(f"Remote path must be absolute: {path}")
        try:
            common = posixpath.commonpath((remote_root, path))
        except ValueError as exc:
            raise ValueError(f"Remote path {path} is outside mapped root {remote_root}") from exc
        if common != remote_root:
            raise ValueError(f"Remote path {path} is outside mapped root {remote_root}")
        relative = posixpath.relpath(path, remote_root)
        if relative == ".":
            return self.local_root
        return self.local_root.joinpath(*PurePosixPath(relative).parts)

    async def map_directory(self, directory):
        return await self.map_path(directory)

    async def _submit(self, request):
        return await self.endpoint._submit(request)

    async def _status(self, directory, run_id):
        return await self.endpoint._status(directory, run_id)

    async def _cancel(self, directory, run_id):
        return await self.endpoint._cancel(directory, run_id)

    async def _follow(self, directory, run_id):
        return await self.endpoint._follow(directory, run_id)

    async def _saved_request(self, directory, run_id):
        return await self.endpoint._saved_request(directory, run_id)

    def create_run(self, reader, spec, callback=None, *, run_id=None):
        return _RemoteRun(self, reader, spec, callback=callback, run_id=run_id)

    async def attach(
        self, reference, reader, callback=None, output_listener=None, start=True, *, request=None
    ):
        expected = f"remote:{self.backend}"
        if reference.executor != expected:
            raise ValueError(
                f"Run reference belongs to {reference.executor!r}, not {expected!r}"
            )
        if not reference.native_id or not reference.directory:
            raise ValueError("Remote run reference lacks run id or directory")

        if request is None:
            status, request = await asyncio.gather(
                self._status(reference.directory, reference.native_id),
                self._saved_request(reference.directory, reference.native_id),
            )
        else:
            status = await self._status(reference.directory, reference.native_id)
        remote_backend = request.get("executor", "local")
        if remote_backend != self.backend:
            raise ValueError(
                f"Remote run uses backend {remote_backend!r}, not {self.backend!r}"
            )

        local_directory = await self.unmap_path(reference.directory)
        local_directory.mkdir(parents=True, exist_ok=True)
        input_remote = status.get("input_file") or reference.input_file
        input_local = await self.unmap_path(input_remote) if input_remote is not None else None

        sync = request.get("client_sync")
        if sync is None:
            sync = {
                "output_file": status.get("output_file") or reference.output_file,
                "error_file": None,
                "output_files": (),
                "cache_files": (),
            }

        output_remote = sync.get("output_file")
        error_remote = sync.get("error_file")
        output_local = await self.unmap_path(output_remote) if output_remote is not None else None
        error_local = await self.unmap_path(error_remote) if error_remote is not None else None

        output_targets = {}
        for remote in sync.get("output_files", ()):
            local = await self.unmap_path(remote)
            output_targets[os.fspath(local)] = remote

        cache_targets = {}
        for item in sync.get("cache_files", ()):
            remote = item["path"]
            before = item.get("before")
            if isinstance(before, list):
                before = tuple(before)
            local = await self.unmap_path(remote)
            cache_targets[os.fspath(local)] = (remote, before)

        output_patterns = tuple(sync.get("output_patterns", ()))
        cache_patterns = tuple(sync.get("cache_patterns", ()))
        cache_before = {}
        for item in sync.get("cache_before", ()):
            before = item.get("signature")
            if isinstance(before, list):
                before = tuple(before)
            cache_before[item["path"]] = before

        opened = []
        try:
            output_file = None
            if output_local is not None:
                output_file = AtomicLocalFile(output_local)
                opened.append(output_file)

            error_file = None
            if error_local is not None:
                if output_local is not None and error_local == output_local:
                    error_file = output_file
                else:
                    error_file = AtomicLocalFile(error_local)
                    opened.append(error_file)

            from .base import ExecutionSpec

            spec = ExecutionSpec(
                argv=(),
                input_file=None if input_local is None else os.fspath(input_local),
                output_file=output_file,
                error_file=error_file,
                directory=os.fspath(local_directory),
                launch=LaunchSpec(
                    local_root=os.fspath(self.local_root),
                    remote_root=self.remote_root,
                    input_sync=self.input_files.value,
                    output_sync=self.output_files.value,
                ),
            )
            return _RemoteRun(
                self,
                reader,
                spec,
                callback=callback,
                output_listener=output_listener,
                start=start,
                reference=reference,
                state=RunState(status.get("state", "unknown")),
                output_targets=output_targets,
                cache_targets=cache_targets,
                output_patterns=output_patterns,
                cache_patterns=cache_patterns,
                cache_before=cache_before,
            )
        except BaseException:
            for file in opened:
                file.discard()
            raise


class RemoteEndpointExecutor(Executor):
    """Reusable remote execution endpoint independent of staging policy."""

    name = "remote"

    def __init__(self, transport, *, backend="local", backend_options=None,
                 command=("ase2sprkkr", "remote")):
        if backend not in _REMOTE_EXECUTORS:
            raise ValueError(f"Unsupported remote executor: {backend!r}")
        self.transport = transport
        self.backend = backend
        self.backend_options = dict(backend_options or {})
        self.command = tuple(command)
        self.computer = getattr(transport, "host", None)
        self._remote_home_cache = None

    async def _remote_home(self):
        if self._remote_home_cache is None:
            result = await self.transport.run(("sh", "-c", 'printf "%s" "$HOME"'))
            home = result.stdout.decode("utf8").strip()
            if not home:
                raise RuntimeError("Could not determine remote home directory")
            self._remote_home_cache = posixpath.normpath(home)
        return self._remote_home_cache

    async def _json_request(self, action, *args, payload=None):
        data = None if payload is None else json.dumps(payload).encode("utf8")
        result = await self.transport.run((*self.command, action, *args), input=data)
        return json.loads(result.stdout.decode("utf8"))

    async def _submit(self, request):
        launch = request.get("launch") or {}
        mpi = launch.get("mpi") or {}
        if (
            self.backend == "local"
            and mpi.get("autodetect")
            and not mpi.get("required")
            and mpi.get("warning", True)
        ):
            resolved = await self._json_request("resolve-mpi", payload={"mpi": "auto"})
            if not resolved.get("mpi"):
                warnings.warn(
                    "MPI autodetection did not find an MPI launcher on the remote host; "
                    "the calculation will run without MPI.",
                    MpiNotFoundWarning,
                    stacklevel=3,
                )
        return await self._json_request("submit", request["directory"], payload=request)

    async def _status(self, directory, run_id):
        return await self._json_request("status", directory, run_id)

    async def _cancel(self, directory, run_id):
        return await self._json_request("cancel", directory, run_id)

    async def _follow(self, directory, run_id):
        return await self.transport.start((*self.command, "follow", directory, run_id))

    async def _saved_request(self, directory, run_id):
        path = posixpath.join(directory, _METADATA_DIR, _RUNS_DIR, run_id, _REQUEST_FILE)
        return json.loads((await self.transport.read_file(path)).decode("utf8"))

    @staticmethod
    def _context(endpoint, launch):
        if not isinstance(launch, LaunchSpec):
            launch = LaunchSpec(**dict(launch))
        return _RemoteRunContext(endpoint, launch)

    def _validate_reference(self, reference):
        expected = f"remote:{self.backend}"
        if reference.executor != expected:
            raise ValueError(
                f"Run reference belongs to {reference.executor!r}, not {expected!r}"
            )
        if not reference.native_id or not reference.directory:
            raise ValueError("Remote run reference lacks run id or directory")

    async def _resolve_submission_directory(self, directory):
        directory = posixpath.normpath(os.fspath(directory))
        if directory == "~":
            return await self._remote_home()
        if directory.startswith("~/"):
            return posixpath.join(await self._remote_home(), directory[2:])
        return directory

    def submission_directory(self, spec):
        launch = spec.launch
        local_root = Path(launch.local_root or Path.home()).expanduser().resolve()
        directory = Path(spec.directory).expanduser().resolve()
        try:
            relative = directory.relative_to(local_root)
        except ValueError as exc:
            raise ValueError(
                f"Local path {directory} is outside mapped root {local_root}; "
                "specify local_root/remote_root explicitly"
            ) from exc
        root = launch.remote_root
        if root is None or root == "~":
            remote_root = "~"
        elif root.startswith("~/"):
            remote_root = posixpath.normpath(root)
        elif posixpath.isabs(root):
            remote_root = posixpath.normpath(root)
        else:
            remote_root = posixpath.normpath(posixpath.join("~", root))
        if not relative.parts:
            return remote_root
        return posixpath.join(remote_root, *PurePosixPath(*relative.parts).parts)

    def status(self, *references: RunReference):
        for reference in references:
            self._validate_reference(reference)
        if not references:
            return ()
        by_directory = {}
        for index, reference in enumerate(references):
            by_directory.setdefault(reference.directory, []).append((index, reference))

        async def query():
            result = [None] * len(references)
            for directory, items in by_directory.items():
                ids = tuple(reference.native_id for _, reference in items)
                states = await self._json_request("status", directory, *ids)
                if len(ids) == 1:
                    states = (states,)
                if len(states) != len(items):
                    raise RuntimeError("Remote endpoint returned wrong number of run statuses")
                for (index, reference), state in zip(items, states):
                    result[index] = RunStatus(
                        reference,
                        RunState(state.get("state", "unknown")),
                        return_code=state.get("return_code"),
                        message=state.get("error"),
                    )
            return tuple(result)

        return _run_sync(query())

    async def _recover_submission(self, run_id, submission):
        directory = await self._resolve_submission_directory(submission.directory)
        state_path = posixpath.join(directory, _METADATA_DIR, _RUNS_DIR, run_id, _STATE_FILE)
        if await self.transport.file_signature(state_path) is None:
            return None
        state = await self._status(directory, run_id)
        reference = RunReference(
            executor=f"remote:{self.backend}", native_id=run_id,
            computer=self.computer, directory=directory,
            input_file=state.get("input_file"), output_file=state.get("output_file"),
            error_file=state.get("error_file"),
        )
        return RunStatus(
            reference=reference, state=RunState(state.get("state", "unknown")),
            return_code=state.get("return_code"), message=state.get("error"),
        )

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        return _run_sync(self._recover_submission(run_id, submission))

    def create_run(self, reader, spec, callback=None):
        return self._context(self, spec.launch).create_run(reader, spec, callback=callback)

    def create_persistent_run(self, reader, spec, *, run_id, callback=None):
        return self._context(self, spec.launch).create_run(
            reader, spec, callback=callback, run_id=run_id
        )

    async def _fetch_result_files(self, reference, files):
        self._validate_reference(reference)
        request = await self._saved_request(reference.directory, reference.native_id)
        launch = LaunchSpec(**(request.get("launch") or {}))
        context = self._context(self, launch)
        remote_root = await context._resolved_remote_root()

        for value in files:
            local = os.path.abspath(os.fspath(value))
            remote = await context.map_path(local)
            if glob.has_magic(local):
                for matched_remote in await self.transport.expand_files(
                    remote, remote_root
                ):
                    matched_local = await context.unmap_path(matched_remote)
                    if not matched_local.is_file():
                        await self.transport.get_file(matched_remote, matched_local)
                continue
            if os.path.isfile(local):
                continue
            if await self.transport.file_signature(remote) is not None:
                await self.transport.get_file(remote, local)

    def fetch_result_files(self, reference, files):
        return _run_sync(self._fetch_result_files(reference, files))

    async def _attach_with_saved_launch(
        self, reference, reader, callback=None, output_listener=None, start=True
    ):
        self._validate_reference(reference)
        request = await self._saved_request(reference.directory, reference.native_id)
        launch = LaunchSpec(**(request.get("launch") or {}))
        context = self._context(self, launch)
        return await context.attach(
            reference, reader, callback=callback, output_listener=output_listener,
            start=start, request=request,
        )

    def attach(
        self, reference: RunReference, reader, callback=None,
        output_listener=None, start=True,
    ):
        return _run_sync(self._attach_with_saved_launch(
            reference, reader, callback=callback,
            output_listener=output_listener, start=start,
        ))
