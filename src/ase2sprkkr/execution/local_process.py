"""Low-level local subprocess launch and stream plumbing."""

from __future__ import annotations

import asyncio
import errno
import os
import signal
import warnings
from dataclasses import replace

from .base import ExecutionSpec, MpiNotFoundWarning
from .mpi import resolve_mpi_command

_DEFAULT_STDIN = object()


def _mpi_command(launch):
    """Materialize portable MPI intent for this machine."""
    mpi = launch.mpi
    command = resolve_mpi_command(mpi)
    if not command:
        if mpi.allows_serial_fallback:
            if mpi.warning:
                warnings.warn(
                    "MPI autodetection did not find an MPI launcher; "
                    "the calculation will run without MPI.",
                    MpiNotFoundWarning,
                    stacklevel=3,
                )
            return None
        raise FileNotFoundError(
            errno.ENOENT,
            "MPI was requested but no MPI launcher was found",
        )
    return command


def materialize_local_spec(spec: ExecutionSpec) -> ExecutionSpec:
    """Choose serial/MPI argv immediately before local process creation."""
    if not spec.launch or spec.mpi_argv is None:
        return spec

    launcher = _mpi_command(spec.launch)
    if launcher is None:
        # ``auto`` may fall back to the serial command on this machine.
        return spec

    return replace(
        spec,
        argv=(*launcher, *spec.mpi_argv),
        stdin=spec.stdin if spec.mpi_uses_stdin else None,
    )


async def start_local_process(
    spec: ExecutionSpec,
    *,
    stdin=_DEFAULT_STDIN,
    stdout=asyncio.subprocess.PIPE,
    stderr=asyncio.subprocess.PIPE,
    start_new_session=False,
):
    """Start a local process with caller-selected standard streams."""
    if stdin is _DEFAULT_STDIN:
        stdin = spec.stdin
    return await asyncio.create_subprocess_exec(
        *spec.argv,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        cwd=spec.directory,
        start_new_session=start_new_session,
    )


async def pump_stream(source, sink, *, eof=None):
    """Forward one subprocess byte stream to a synchronous chunk consumer."""
    try:
        while data := await source.read(65536):
            sink(data)
    finally:
        if eof is not None:
            eof()


async def stop_local_process(
    process, *, process_group=False, timeout=0.1, waiting=None,
):
    """Send SIGTERM, then SIGKILL on timeout, without cancelling process.wait().

    For a process group, signal even if its leader has already exited:
    descendants may still be alive and need the full SIGTERM grace period.
    ``waiting`` reuses an existing wait task.
    """
    if process is None or (process.returncode is not None and not process_group):
        return

    if waiting is None:
        waiting = asyncio.create_task(process.wait())

    try:
        if process_group:
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        await waiting
        return

    deadline = asyncio.get_running_loop().time() + timeout
    try:
        await asyncio.wait_for(asyncio.shield(waiting), timeout)
    except asyncio.TimeoutError:
        if process.returncode is not None and not process_group:
            return
    else:
        if not process_group:
            return
        # The group can outlive its leader. Preserve the SIGTERM grace
        # period if any descendant still belongs to the group.
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return
        await asyncio.sleep(max(0, deadline - asyncio.get_running_loop().time()))

    try:
        if process_group:
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
    except ProcessLookupError:
        pass
    await waiting
