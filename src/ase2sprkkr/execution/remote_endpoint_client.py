"""Client-owned remote endpoint lifecycle.

This module keeps staging and file synchronization from :mod:`remote_endpoint`,
but removes the target-side ASE2SPRKKR run database.  The target returns one
opaque ``remote_id``; the client persists that id as ``RunReference.native_id``
and keeps reconnect-only staging information in ``RunReference.context``.
"""

from __future__ import annotations

import glob
import os
from uuid import uuid4

from .base import RunReference, RunState, RunStatus, SubmissionNotStarted
from .remote import _run_sync
from .remote_endpoint import (
    RemoteEndpointExecutor as _LegacyRemoteEndpointExecutor,
    _RemoteRun as _LegacyRemoteRun,
    _RemoteRunContext as _LegacyRemoteRunContext,
)


def _reference_context(reference):
    context = reference.context or {}
    if not isinstance(context, dict):
        raise ValueError("Remote run reference has invalid reconnect context")
    return context


class _ClientOwnedRemoteRun(_LegacyRemoteRun):
    """Remote run whose reconnect data is persisted only by the client."""

    async def _submit_async(self):
        run_id = self._run_id or uuid4().hex
        try:
            remote_directory = await self.executor.map_directory(self.spec.directory)
            await self.executor.transport.run(("mkdir", "-p", "--", remote_directory))
            request = await self._prepare_request(run_id, remote_directory)
        except Exception as exc:
            raise SubmissionNotStarted(str(exc)) from exc

        status = await self.executor._submit(request)
        remote_id = status.get("remote_id")
        if not remote_id:
            raise RuntimeError("Remote submit returned no remote_id")
        context = {
            "launch": request.get("launch") or {},
            "client_sync": request.get("client_sync") or {},
        }
        return remote_id, context, status

    def submit(self):
        with self._lock:
            if self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self._state.value}")

            remote_id, context, status = _run_sync(self._submit_async())
            self._state = RunState(status["state"])
            self._return_code = status.get("return_code")
            self._reference = RunReference(
                executor=f"remote:{self.executor.backend}",
                native_id=str(remote_id),
                computer=self.executor.computer,
                directory=status["directory"],
                input_file=status.get("input_file"),
                output_file=status.get("output_file"),
                error_file=status.get("error_file"),
                context=context,
            )
        self._start_monitor()


class _ClientOwnedRunContext(_LegacyRemoteRunContext):
    def create_run(self, reader, spec, callback=None, *, run_id=None):
        return _ClientOwnedRemoteRun(
            self,
            reader,
            spec,
            callback=callback,
            run_id=run_id,
        )

    async def attach(
        self,
        reference,
        reader,
        callback=None,
        output_listener=None,
        start=True,
        *,
        request=None,
    ):
        context = _reference_context(reference)
        request = {
            "executor": self.backend,
            "launch": context.get("launch") or {},
            "client_sync": context.get("client_sync") or {},
        }
        return await super().attach(
            reference,
            reader,
            callback=callback,
            output_listener=output_listener,
            start=start,
            request=request,
        )


class RemoteEndpointExecutor(_LegacyRemoteEndpointExecutor):
    """Remote endpoint with one client-side lifecycle authority.

    The target never owns an ASE2SPRKKR ``RunState``.  It only accepts opaque
    backend ids for submit/status/cancel/follow operations.
    """

    @staticmethod
    def _context(endpoint, launch):
        from .base import LaunchSpec

        if not isinstance(launch, LaunchSpec):
            launch = LaunchSpec(**dict(launch))
        return _ClientOwnedRunContext(endpoint, launch)

    async def _submit(self, request):
        launch = request.get("launch") or {}
        mpi = launch.get("mpi") or {}
        if (
            self.backend == "local"
            and mpi.get("autodetect")
            and not mpi.get("required")
            and mpi.get("warning", True)
        ):
            from .base import MpiNotFoundWarning
            import warnings

            resolved = await self._json_request("resolve-mpi", payload={"mpi": "auto"})
            if not resolved.get("mpi"):
                warnings.warn(
                    "MPI autodetection did not find an MPI launcher on the remote host; "
                    "the calculation will run without MPI.",
                    MpiNotFoundWarning,
                    stacklevel=3,
                )
        return await self._json_request("submit", request["directory"], payload=request)

    async def _status(self, directory, native_id):
        return await self._json_request("status", native_id)

    async def _cancel(self, directory, native_id):
        return await self._json_request("cancel", native_id)

    async def _follow(self, directory, native_id):
        from .stream_protocol import FramedRunningCommand

        process = await self.transport.start(
            (*self.command, "follow", str(native_id), "--framed")
        )
        return FramedRunningCommand(process)

    async def _saved_request(self, directory, run_id):
        raise RuntimeError("Remote lifecycle metadata is client-owned")

    def status(self, *references: RunReference):
        for reference in references:
            self._validate_reference(reference)
        if not references:
            return ()

        async def query():
            result = []
            for reference in references:
                state = await self._status(reference.directory, reference.native_id)
                result.append(
                    RunStatus(
                        reference,
                        RunState(state.get("state", "unknown")),
                        return_code=state.get("return_code"),
                        message=state.get("error"),
                    )
                )
            return tuple(result)

        return _run_sync(query())

    async def _recover_submission(self, run_id, submission, *, submitted_at=None):
        directory = await self._resolve_submission_directory(submission.directory)
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
            raise NotImplementedError(state.get("error") or "Remote recovery unsupported")
        remote_id = state.get("remote_id")
        if not remote_id:
            raise RuntimeError("Recovered remote submission has no remote_id")

        reference = RunReference(
            executor=f"remote:{self.backend}",
            native_id=str(remote_id),
            computer=self.computer,
            directory=state.get("directory") or directory,
            input_file=state.get("input_file"),
            output_file=state.get("output_file"),
            error_file=state.get("error_file"),
            context=getattr(submission, "context", None),
        )
        return RunStatus(
            reference,
            RunState(state.get("state", "unknown")),
            return_code=state.get("return_code"),
            message=state.get("error"),
        )

    def recover_submission(self, run_id, submission, *, submitted_at=None):
        return _run_sync(
            self._recover_submission(run_id, submission, submitted_at=submitted_at)
        )

    async def _fetch_result_files(self, reference, files):
        self._validate_reference(reference)
        context_data = _reference_context(reference)
        from .base import LaunchSpec

        launch = LaunchSpec(**(context_data.get("launch") or {}))
        context = self._context(self, launch)
        remote_root = await context._resolved_remote_root()

        for value in files:
            local = os.path.abspath(os.fspath(value))
            remote = await context.map_path(local)
            if glob.has_magic(local):
                for matched_remote in await self.transport.expand_files(remote, remote_root):
                    matched_local = await context.unmap_path(matched_remote)
                    if not matched_local.is_file():
                        await self.transport.get_file(matched_remote, matched_local)
                continue
            if os.path.isfile(local):
                continue
            if await self.transport.file_signature(remote) is not None:
                await self.transport.get_file(remote, local)

    async def _attach_with_saved_launch(
        self,
        reference,
        reader,
        callback=None,
        output_listener=None,
        start=True,
    ):
        self._validate_reference(reference)
        context_data = _reference_context(reference)
        from .base import LaunchSpec

        launch = LaunchSpec(**(context_data.get("launch") or {}))
        context = self._context(self, launch)
        return await context.attach(
            reference,
            reader,
            callback=callback,
            output_listener=output_listener,
            start=start,
        )
