"""Endpoint-only remote executor with per-run staging configuration."""

from __future__ import annotations

from dataclasses import asdict

from .base import LaunchSpec, RunReference
from .remote import RemoteExecutor, _run_sync


class _RemoteRunContext(RemoteExecutor):
    """Lightweight configured facade for one remote run.

    It shares transport/backend identity with the pooled endpoint executor but
    keeps filesystem mapping and sync policy local to the concrete run.
    """

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
        # Reuse the endpoint's cached remote HOME lookup when available.
        self._remote_home_cache = endpoint._remote_home_cache

    async def _remote_home(self):
        home = await self.endpoint._remote_home()
        self._remote_home_cache = home
        return home


class RemoteEndpointExecutor(RemoteExecutor):
    """Reusable remote endpoint independent of launch/staging policy."""

    def __init__(
        self,
        transport,
        *,
        backend="local",
        backend_options=None,
        command=("ase2sprkkr", "remote"),
    ):
        super().__init__(
            transport,
            backend=backend,
            backend_options=backend_options,
            command=command,
        )

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

    def launch_snapshot(self, spec):
        """Return the serializable launch configuration carried by *spec*."""
        return asdict(spec.launch)
