"""Backend status queries for persistent execution references."""

from __future__ import annotations

import errno
import os
from functools import singledispatch

from .base import Executor, RunReference, RunState, RunStatus
from .local import LocalExecutor
from .psij import PsijExecutor, _STATE_MAP
from .remote import RemoteExecutor, _run_sync
from .remote_endpoint import RemoteEndpointExecutor


def _validate_executor(references, expected):
    for reference in references:
        if reference.executor != expected:
            raise ValueError(
                "Run reference belongs to {!r}, not {!r}".format(
                    reference.executor, expected
                )
            )


@singledispatch
def status(executor: Executor, *references: RunReference):
    """Return ordered status snapshots for references handled by *executor*."""
    if not references:
        return ()
    raise NotImplementedError(
        "{} cannot query persistent run status".format(type(executor).__name__)
    )


def _pid_running(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        if exc.errno == errno.ESRCH:
            return False
        if exc.errno == errno.EPERM:
            return True
        raise
    return True


@status.register(LocalExecutor)
def _(executor, *references):
    _validate_executor(references, executor.name)
    out = []
    for reference in references:
        try:
            pid = int(reference.native_id)
        except (TypeError, ValueError):
            state = RunState.UNKNOWN
        else:
            # Local execution has no reconnectable controller. Once the recorded
            # process is gone, its final outcome cannot be inferred from the PID.
            state = RunState.RUNNING if _pid_running(pid) else RunState.UNKNOWN
        out.append(RunStatus(reference, state))
    return tuple(out)


@status.register(PsijExecutor)
def _(executor, *references):
    expected = "psij:{}".format(executor.backend)
    _validate_executor(references, expected)
    jobs = getattr(executor, "_status_jobs", None)
    if jobs is None:
        jobs = executor._status_jobs = {}

    out = []
    for reference in references:
        attached = False
        job = jobs.get(reference.native_id)
        if job is None:
            job = executor._psij.Job()
            executor._job_executor.attach(job, reference.native_id)
            jobs[reference.native_id] = job
            attached = True
        status_value = getattr(job, "status", None)
        state_value = getattr(status_value, "state", None)
        state_name = getattr(state_value, "name", None)
        state = (
            RunState.UNKNOWN
            if attached and state_name == "NEW"
            else _STATE_MAP.get(state_name, RunState.UNKNOWN)
        )
        out.append(
            RunStatus(
                reference,
                state,
                return_code=getattr(status_value, "exit_code", None),
                message=getattr(status_value, "message", None),
            )
        )
    return tuple(out)


@status.register(RemoteEndpointExecutor)
@status.register(RemoteExecutor)
def _(executor, *references):
    expected = "remote:{}".format(executor.backend)
    _validate_executor(references, expected)
    if not references:
        return ()

    by_directory = {}
    for index, reference in enumerate(references):
        if not reference.directory:
            raise ValueError("Remote run reference lacks directory")
        by_directory.setdefault(reference.directory, []).append((index, reference))

    async def query():
        result = [None] * len(references)
        for directory, items in by_directory.items():
            ids = tuple(reference.native_id for _, reference in items)
            states = await executor._json_request("status", directory, *ids)
            if len(ids) == 1:
                states = (states,)
            for (index, reference), state_value in zip(items, states):
                result[index] = RunStatus(
                    reference,
                    RunState(state_value.get("state", "unknown")),
                    return_code=state_value.get("return_code"),
                    message=state_value.get("error"),
                )
        return tuple(result)

    return _run_sync(query())
