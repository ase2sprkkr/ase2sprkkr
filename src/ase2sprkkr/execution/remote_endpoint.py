"""Reusable remote endpoint with per-run staging configuration."""

from __future__ import annotations

import json
import posixpath

from .base import Executor, LaunchSpec, RunReference
from .remote import (
    RemoteExecutor as _RemoteStagingExecutor,
    _METADATA_DIR,
    _REQUEST_FILE,
    _REMOTE_EXECUTORS,
    _RUNS_DIR,
    _run_sync,
)


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

    def create_run(self, reader, spec, callback=None):
        """Create one run using staging parameters embedded in *spec*."""
        return self._context(self, spec.launch).create_run(
            reader, spec, callback=callback
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
