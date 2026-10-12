"""Reusable remote endpoint with per-run staging configuration."""

from __future__ import annotations

import json
import os
import posixpath
import warnings
from pathlib import Path, PurePosixPath
from uuid import uuid4

from .base import Executor, LaunchSpec, MpiNotFoundWarning, RunReference
from .remote import (
    RemoteExecutor as _RemoteStagingExecutor,
    RemoteRunner,
    _METADATA_DIR,
    _REQUEST_FILE,
    _REMOTE_EXECUTORS,
    _RUNS_DIR,
    _run_sync,
)


class _RemoteRun(RemoteRunner):
    """Remote runner that can reuse an identity chosen before submission."""

    def __init__(self, *args, run_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._run_id = run_id

    async def _submit_async(self):
        run_id = self._run_id or uuid4().hex
        remote_directory = await self.executor.map_directory(self.spec.directory)
        await self.executor.transport.run(
            ("mkdir", "-p", "--", remote_directory)
        )
        request = await self._prepare_request(run_id, remote_directory)
        status = await self.executor._submit(request)
        return run_id, status


class _RemoteRunContext(_RemoteStagingExecutor):
    """Per-run staging facade sharing one pooled endpoint transport."""

    def __init__(self, endpoint, launch):
        self.endpoint = endpoint
        super().__init__(
            endpoint.transport,
            backend=endpoint.backend,
            backend_options=endpoint.backend_options,
            local_root=launch.local_root,
            remote_root=launch.remote_root,
            command=endpoint.command,
            input_files=launch.input_sync,
            output_files=launch.output_sync,
        )

    async def _remote_home(self):
        return await self.endpoint._remote_home()

    async def _submit(self, request):
        return await self.endpoint._submit(request)

    def create_run(self, reader, spec, callback=None, *, run_id=None):
        return _RemoteRun(
            self,
            reader,
            spec,
            callback=callback,
            run_id=run_id,
        )


class RemoteEndpointExecutor(Executor):
    """Reusable remote execution endpoint independent of staging policy."""

    name = "remote"

    def __init__(
        self,
        transport,
        *,
        backend="local",
        backend_options=None,
        command=("ase2sprkkr", "remote"),
    ):
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
            result = await self.transport.run(
                ("sh", "-c", 'printf "%s" "$HOME"')
            )
            home = result.stdout.decode("utf8").strip()
            if not home:
                raise RuntimeError("Could not determine remote home directory")
            self._remote_home_cache = posixpath.normpath(home)
        return self._remote_home_cache

    async def _json_request(self, action, *args, payload=None):
        data = None if payload is None else json.dumps(payload).encode("utf8")
        result = await self.transport.run(
            (*self.command, action, *args), input=data
        )
        return json.loads(result.stdout.decode("utf8"))

    async def _submit(self, request):
        launch = request.get("launch") or {}
        if (
            self.backend == "local"
            and launch.get("mpi") == "auto"
            and launch.get("mpi_warning", True)
        ):
            resolved = await self._json_request(
                "resolve-mpi",
                payload={"mpi": "auto"},
            )
            if not resolved.get("mpi"):
                warnings.warn(
                    "MPI autodetection did not find an MPI launcher on the remote host; "
                    "the calculation will run without MPI.",
                    MpiNotFoundWarning,
                    stacklevel=3,
                )
        return await self._json_request(
            "submit", request["directory"], payload=request
        )

    async def _status(self, directory, run_id):
        return await self._json_request("status", directory, run_id)

    async def _cancel(self, directory, run_id):
        return await self._json_request("cancel", directory, run_id)

    async def _follow(self, directory, run_id):
        return await self.transport.start(
            (*self.command, "follow", directory, run_id)
        )

    async def _saved_request(self, directory, run_id):
        path = posixpath.join(
            directory,
            _METADATA_DIR,
            _RUNS_DIR,
            run_id,
            _REQUEST_FILE,
        )
        return json.loads((await self.transport.read_file(path)).decode("utf8"))

    @staticmethod
    def _context(endpoint, launch):
        if not isinstance(launch, LaunchSpec):
            launch = LaunchSpec(**dict(launch))
        return _RemoteRunContext(endpoint, launch)

    def submission_directory(self, spec):
        """Return a stable target-side namespace without contacting the host."""
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

    def create_run(self, reader, spec, callback=None):
        """Create a nonpersistent run with a backend-generated identity."""
        return self._context(self, spec.launch).create_run(
            reader, spec, callback=callback
        )

    def create_persistent_run(self, reader, spec, *, run_id, callback=None):
        """Create a remote run using the manager's persistent run identity."""
        return self._context(self, spec.launch).create_run(
            reader, spec, callback=callback, run_id=run_id
        )

    async def _attach_with_saved_launch(self, reference, reader, callback=None):
        request = await self._saved_request(reference.directory, reference.native_id)
        launch = LaunchSpec(**(request.get("launch") or {}))
        context = self._context(self, launch)
        return await context._attach(reference, reader, callback=callback)

    def attach(self, reference: RunReference, reader, callback=None):
        """Reconnect using the run's saved launch/staging parameters."""
        return _run_sync(
            self._attach_with_saved_launch(reference, reader, callback=callback)
        )
