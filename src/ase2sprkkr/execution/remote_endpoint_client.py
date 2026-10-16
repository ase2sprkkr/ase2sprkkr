"""Client-owned remote endpoint lifecycle.

This module keeps staging and file synchronization from :mod:`remote_endpoint`,
but removes the target-side ASE2SPRKKR run database.  The target returns one
opaque ``remote_id``; the client persists that id together with the per-run
reconnect context in :class:`RunReference.native_id`.
"""

from __future__ import annotations

import base64
import glob
import json
import os
from dataclasses import asdict
from uuid import uuid4

from .base import RunReference, RunState, RunStatus, SubmissionNotStarted
from .remote import _run_sync
from .remote_endpoint import (
    RemoteEndpointExecutor as _LegacyRemoteEndpointExecutor,
    _RemoteRun as _LegacyRemoteRun,
    _RemoteRunContext as _LegacyRemoteRunContext,
)


def _pack_reference(remote_id, *, launch, client_sync):
    payload = {
        "remote_id": str(remote_id),
        "launch": launch,
        "client_sync": client_sync,
    }
    raw = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf8")
    return "v1:" + base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def _unpack_reference(native_id):
    value = str(native_id)
    if not value.startswith("v1:"):
        raise ValueError("Remote run reference does not contain client reconnect data")
    value = value[3:]
    padding = "=" * (-len(value) % 4)
    try:
        payload = json.loads(base64.urlsafe_b64decode(value + padding).decode("utf8"))
    except (ValueError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid remote run reference") from exc
    if not isinstance(payload, dict) or not payload.get("remote_id"):
        raise ValueError("Invalid remote run reference")
    return payload


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
        native_id = _pack_reference(
            remote_id,
            launch=request.get("launch") or {},
            client_sync=request.get("client_sync") or {},
        )
        return native_id, status


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
        packed = _unpack_reference(reference.native_id)
        request = {
            "executor": self.backend,
            "launch": packed.get("launch") or {},
            "client_sync": packed.get("client_sync") or {},
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

    @staticmethod
    def _remote_id(reference):
        return _unpack_reference(reference.native_id)["remote_id"]

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
        remote_id = _unpack_reference(native_id)["remote_id"]
        return await self._json_request("status", remote_id)

    async def _cancel(self, directory, native_id):
        remote_id = _unpack_reference(native_id)["remote_id"]
        return await self._json_request("cancel", remote_id)

    async def _follow(self, directory, native_id):
        from .stream_protocol import FramedRunningCommand

        remote_id = _unpack_reference(native_id)["remote_id"]
        process = await self.transport.start((*self.command, "follow", remote_id, "--framed"))
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
        try:
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
        except RuntimeError:
            raise
        if state is None:
            return None
        if state.get("unsupported"):
            raise NotImplementedError(state.get("error") or "Remote recovery unsupported")
        remote_id = state.get("remote_id")
        if not remote_id:
            raise RuntimeError("Recovered remote submission has no remote_id")

        # Recovery happens before a RunReference exists.  The persisted client
        # record still contains result manifests, while launch/file-sync context
        # is unavailable unless the backend was able to return a reference in
        # the original submit response.  Scheduler recovery is therefore kept
        # deliberately minimal here; normal attach uses the reference created
        # by a successful submit.
        native_id = _pack_reference(remote_id, launch={}, client_sync={})
        reference = RunReference(
            executor=f"remote:{self.backend}",
            native_id=native_id,
            computer=self.computer,
            directory=state.get("directory") or directory,
            input_file=state.get("input_file"),
            output_file=state.get("output_file"),
            error_file=state.get("error_file"),
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
        packed = _unpack_reference(reference.native_id)
        from .base import LaunchSpec

        launch = LaunchSpec(**(packed.get("launch") or {}))
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
        packed = _unpack_reference(reference.native_id)
        from .base import LaunchSpec

        launch = LaunchSpec(**(packed.get("launch") or {}))
        context = self._context(self, launch)
        return await context.attach(
            reference,
            reader,
            callback=callback,
            output_listener=output_listener,
            start=start,
        )
