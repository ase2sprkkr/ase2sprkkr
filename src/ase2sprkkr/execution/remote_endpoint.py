"""Remote execution with client-owned lifecycle persistence."""

from __future__ import annotations

import glob
import json
import os
import posixpath
import warnings
from dataclasses import asdict, replace
from pathlib import Path, PurePosixPath
from uuid import uuid4

from .base import (
    ExecutionSpec,
    Executor,
    FailurePhase,
    LaunchSpec,
    MpiNotFoundWarning,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .remote import FileSyncPolicy, RemoteRunner, _run_sync
from .stream_protocol import FramedRunningCommand
from .transport import AtomicLocalFile
from .run_model import FileEntry, CacheEntry, FileManifest, FileMapping, WorkingDirectory
from .endpoint import endpoint_for_executor


_REMOTE_BACKENDS = {"local", "psij"}


def _local_path(value, directory):
    if value is None:
        return None
    value = getattr(value, "name", value)
    try:
        path = os.fspath(value)
    except TypeError:
        return None
    if isinstance(path, int):
        return None
    if not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


def _client_status_state(state):
    """Translate target backend facts to the remote executor's client lifecycle."""
    backend_state = RunState(state.get("state", "unknown"))
    raw_phase = state.get("failure_phase")
    phase = None if raw_phase is None else FailurePhase(raw_phase)

    if backend_state == RunState.FINISHED:
        return RunState.DOWNLOADING
    if backend_state == RunState.FAILED and phase not in (
        FailurePhase.SUBMISSION,
        FailurePhase.LAUNCH,
    ):
        return RunState.DOWNLOADING
    return backend_state


class _RemoteRun(RemoteRunner):
    """One remote run whose persistent lifecycle is owned by the client."""

    def __init__(self, *args, run_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._run_id = run_id

    async def _monitor(self):
        local_output = self.result.output_file
        result = await super()._monitor()
        remote_output = None if self.files is None or self.files.stdout is None else self.files.stdout.remote_path
        if (
            self.state == RunState.FINISHED
            and local_output is not None
            and remote_output is not None
        ):
            signature = await self.executor.transport.file_signature(remote_output)
            if signature is not None:
                stat = os.stat(local_output)
                if stat.st_size != signature[0]:
                    self._state = RunState.FAILED
                    self._failure_phase = FailurePhase.DOWNLOAD
                    self._failure_message = (
                        "Streamed local output size does not match the authoritative "
                        f"remote output: {local_output}"
                    )
                    raise RuntimeError(self._failure_message)
                os.utime(local_output, ns=(stat.st_atime_ns, signature[1]))
        return result

    async def _submitted_manifest(self, request):
        directory = os.path.abspath(os.fspath(self.spec.directory))
        def local(value):
            return _local_path(value, directory)
        inputs = []
        for item in (self.spec.input_file, *self.spec.input_files):
            path = local(item)
            if path is not None:
                inputs.append(FileEntry(path, await self.executor.map_path(path)))
        outputs = []
        for item in self.spec.output_files:
            path = local(item)
            if path is not None and path not in (self.result.output_file, self.result.error_file):
                outputs.append(FileEntry(path, await self.executor.map_path(path)))
        caches = []
        for item in self.spec.cache_files:
            path = local(item)
            if path is None:
                continue
            remote = await self.executor.map_path(path)
            previous = self._cache_targets.get(path)
            caches.append(CacheEntry(FileEntry(path, remote), None if previous is None else previous[1]))
        for remote, before in self._cache_before.items():
            caches.append(CacheEntry(FileEntry(remote_path=remote), before))
        sync = request.get("client_sync") or {}
        mapping = FileMapping(
            local_root=os.fspath(self.executor.local_root),
            remote_root=await self.executor._resolved_remote_root(),
            input_sync=self.executor.input_files.value,
            output_sync=self.executor.output_files.value,
        )
        return FileManifest(
            inputs=tuple(inputs), outputs=tuple(outputs), caches=tuple(caches),
            stdout=None if self.result.output_file is None else FileEntry(self.result.output_file, sync.get("output_file")),
            stderr=None if self.result.error_file is None else FileEntry(self.result.error_file, sync.get("error_file")),
            mapping=mapping,
            availability_state="pending",
        )

    async def _submit_async(self):
        run_id = self._run_id or uuid4().hex
        try:
            remote_directory = await self.executor.map_directory(self.spec.directory)
            await self.executor.transport.run(("mkdir", "-p", "--", remote_directory))
            request = await self._prepare_request(run_id, remote_directory)
        except Exception as exc:
            self._failure_phase = FailurePhase.SUBMISSION
            self._failure_message = str(exc)
            raise SubmissionNotStarted(str(exc)) from exc

        status = await self.executor._submit(request)
        remote_id = status.get("remote_id")
        if not remote_id:
            raise RuntimeError("Remote submit returned no remote_id")
        manifest = await self._submitted_manifest(request)
        return remote_id, manifest, status, request

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")
            remote_id, manifest, status, request = _run_sync(self._submit_async())
            self._state = RunState(status["state"])
            self._return_code = status.get("return_code")
            self.result.return_code = self._return_code
            phase = status.get("failure_phase")
            self._failure_phase = None if phase is None else FailurePhase(phase)
            self._failure_message = status.get("error")
            self._backend_id = str(remote_id)
            self._remote_directory = request["directory"]
            self._files = manifest
        self._start_monitor()


class _RemoteRunContext:
    """Per-run filesystem and synchronization policy."""

    def __init__(self, endpoint, launch):
        self.endpoint = endpoint
        self.transport = endpoint.transport
        self.backend = endpoint.backend
        self.backend_options = endpoint.backend_options
        self.command = endpoint.command
        self.computer = endpoint.computer
        self.local_root = Path(launch.local_root or Path.home()).expanduser().resolve()
        self.remote_root = (
            None if launch.remote_root is None else os.fspath(launch.remote_root)
        )
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
            raise ValueError(
                f"Remote path {path} is outside mapped root {remote_root}"
            ) from exc
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

    async def _status(self, directory, backend_id, files=None):
        return await self.endpoint._status(directory, backend_id, files)

    async def _cancel(self, directory, backend_id, files=None):
        return await self.endpoint._cancel(directory, backend_id, files)

    async def _follow(self, directory, backend_id, files=None):
        return await self.endpoint._follow(directory, backend_id, files)

    def create_run(self, spec, callback=None, *, run_id=None):
        return _RemoteRun(self, spec, callback=callback, run_id=run_id)

    @staticmethod
    def _manifest_sync(files):
        return {
            "output_files": [
                entry.remote_path for entry in files.outputs
                if entry.remote_path and not glob.has_magic(entry.remote_path)
            ],
            "output_patterns": [
                entry.remote_path for entry in files.outputs
                if entry.remote_path and glob.has_magic(entry.remote_path)
            ],
            "cache_files": [
                {"path": item.file.remote_path, "before": item.before}
                for item in files.caches
                if item.file.local_path is not None and item.file.remote_path
                and not glob.has_magic(item.file.remote_path)
            ],
            "cache_patterns": [
                item.file.remote_path for item in files.caches
                if item.file.local_path is not None and item.file.remote_path
                and glob.has_magic(item.file.remote_path)
            ],
            "cache_before": [
                {"path": item.file.remote_path, "signature": item.before}
                for item in files.caches if item.file.local_path is None
            ],
        }

    async def attach(self, run, callback=None, output_listener=None, start=True):
        execution = run.execution
        if execution.backend_id is None or execution.directory.remote is None:
            raise ValueError("Remote run lacks backend ID or target directory")
        status = await self._status(execution.directory.remote, execution.backend_id, run.files)
        sync = self._manifest_sync(run.files)
        local_directory = Path(execution.directory.local)
        local_directory.mkdir(parents=True, exist_ok=True)
        input_local = run.files.input_file
        output_local = run.files.output_file
        error_local = None if run.files.stderr is None else run.files.stderr.local_path

        output_targets = {}
        for remote in sync.get("output_files", ()):
            local = await self.unmap_path(remote)
            output_targets[os.fspath(local)] = remote

        cache_targets = {}
        for item in sync.get("cache_files", ()):
            before = item.get("before")
            if isinstance(before, list):
                before = tuple(before)
            local = await self.unmap_path(item["path"])
            cache_targets[os.fspath(local)] = (item["path"], before)

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
            phase = status.get("failure_phase")
            runner = _RemoteRun(
                self,
                spec,
                callback=callback,
                output_listener=output_listener,
                start=False,
                backend_id=execution.backend_id,
                remote_directory=execution.directory.remote,
                files=run.files,
                state=_client_status_state(status),
                output_targets=output_targets,
                cache_targets=cache_targets,
                output_patterns=tuple(sync.get("output_patterns", ())),
                cache_patterns=tuple(sync.get("cache_patterns", ())),
                cache_before=cache_before,
            )
            runner._return_code = status.get("return_code")
            runner.result.return_code = runner._return_code
            runner._failure_phase = None if phase is None else FailurePhase(phase)
            runner._failure_message = status.get("error")
            if start:
                runner.start()
            return runner
        except BaseException:
            for file in opened:
                file.discard()
            raise


class RemoteEndpointExecutor(Executor):
    """Reusable remote endpoint with client-owned lifecycle state."""

    name = "remote"

    def __init__(
        self,
        transport,
        *,
        backend="local",
        backend_options=None,
        command=("ase2sprkkr", "remote"),
    ):
        if backend not in _REMOTE_BACKENDS:
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

    @staticmethod
    def _context(endpoint, launch):
        if not isinstance(launch, LaunchSpec):
            launch = LaunchSpec(**dict(launch))
        return _RemoteRunContext(endpoint, launch)

    def _validate_run(self, run):
        if run.execution.endpoint != endpoint_for_executor(self):
            raise ValueError("Run belongs to a different remote endpoint")
        if not run.execution.backend_id or not run.execution.directory.remote:
            raise ValueError("Remote run lacks backend identity or working directory")

    @staticmethod
    def _launch_from_manifest(files):
        mapping = files.mapping
        if mapping is None:
            raise ValueError("Remote run has no persisted file mapping")
        return LaunchSpec(
            local_root=mapping.local_root, remote_root=mapping.remote_root,
            input_sync=mapping.input_sync, output_sync=mapping.output_sync,
        )

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

    async def _submit(self, request):
        launch = request.get("launch") or {}
        mpi = launch.get("mpi") or {}
        if (
            self.backend == "local"
            and mpi.get("autodetect")
            and not mpi.get("required")
            and mpi.get("warning", True)
        ):
            resolved = await self._json_request(
                "resolve-mpi", payload={"mpi": "auto"}
            )
            if not resolved.get("mpi"):
                warnings.warn(
                    "MPI autodetection did not find an MPI launcher on the remote "
                    "host; the calculation will run without MPI.",
                    MpiNotFoundWarning,
                    stacklevel=3,
                )
        return await self._json_request(
            "submit", request["directory"], payload=request
        )

    async def _status(self, directory, backend_id, files=None):
        return await self._json_request("status", backend_id)

    async def _cancel(self, directory, backend_id, files=None):
        return await self._json_request("cancel", backend_id)

    async def _follow(self, directory, backend_id, files=None):
        options = []
        if files is not None:
            if files.stdout and files.stdout.remote_path:
                options.extend(("--stdout-path", files.stdout.remote_path))
            if files.stderr and files.stderr.remote_path:
                options.extend(("--stderr-path", files.stderr.remote_path))
        process = await self.transport.start(
            (*self.command, "follow", str(backend_id), *options, "--framed")
        )
        return FramedRunningCommand(process)

    def acknowledge(self, run):
        """Release a delivered receipt only after durable RunRecord completion."""
        self._validate_run(run)
        return _run_sync(self._json_request("ack", run.execution.backend_id))

    def status(self, *runs):
        for run in runs:
            self._validate_run(run)
        if not runs:
            return ()

        async def query():
            result = []
            for run in runs:
                execution = run.execution
                state = await self._status(
                    execution.directory.remote, execution.backend_id, run.files
                )
                result.append(RunStatus(
                    _client_status_state(state),
                    backend_id=execution.backend_id,
                    return_code=state.get("return_code"),
                    message=state.get("error"),
                    failure_phase=state.get("failure_phase"),
                ))
            return tuple(result)

        return _run_sync(query())

    async def _recover_submission(self, run_id, submission, *, submitted_at=None):
        directory = await self._resolve_submission_directory(submission.execution.directory.target)
        state = await self._json_request(
            "recover",
            directory,
            run_id,
            payload={
                "executor": self.backend,
                "executor_options": self.backend_options,
                "submitted_at": submitted_at,
            },
        )
        if state is None:
            return None
        if state.get("unsupported"):
            raise NotImplementedError(
                state.get("error") or "Remote recovery unsupported"
            )
        remote_id = state.get("remote_id")
        if not remote_id:
            raise RuntimeError("Recovered remote submission has no remote_id")

        execution = submission.execution.submitted(remote_id)
        target_directory = state.get("directory") or directory
        execution = replace(
            execution, directory=WorkingDirectory(
                execution.directory.local, target_directory,
            ),
        )
        files = submission.files
        changes = {}
        for field, key in (("stdout", "output_file"), ("stderr", "error_file")):
            remote = state.get(key)
            if remote is not None:
                old = getattr(files, field)
                changes[field] = FileEntry(
                    None if old is None else old.local_path, remote,
                )
        if changes:
            files = replace(files, **changes)
        return RunStatus(
            _client_status_state(state),
            backend_id=str(remote_id),
            execution=execution,
            files=files,
            return_code=state.get("return_code"),
            message=state.get("error"),
            failure_phase=state.get("failure_phase"),
        )

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        return _run_sync(
            self._recover_submission(
                run_id, submission, submitted_at=submitted_at
            )
        )

    def create_run(self, spec, callback=None):
        return self._context(self, spec.launch).create_run(spec, callback=callback)

    def create_persistent_run(self, spec, *, run_id, callback=None):
        return self._context(self, spec.launch).create_run(
            spec, callback=callback, run_id=run_id
        )

    async def _fetch_result_files(self, reference, files):
        self._validate_reference(reference)
        context_data = _reference_context(reference)
        launch = LaunchSpec(**(context_data.get("launch") or {}))
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

    async def _attach(
        self,
        reference,
        callback=None,
        output_listener=None,
        start=True,
    ):
        self._validate_reference(reference)
        context_data = _reference_context(reference)
        launch = LaunchSpec(**(context_data.get("launch") or {}))
        return await self._context(self, launch).attach(
            reference,
            callback=callback,
            output_listener=output_listener,
            start=start,
        )

    def attach(
        self,
        reference: RunReference,
        callback=None,
        output_listener=None,
        start=True,
    ):
        return _run_sync(
            self._attach(
                reference,
                callback=callback,
                output_listener=output_listener,
                start=start,
            )
        )
