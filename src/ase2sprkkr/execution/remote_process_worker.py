"""Direct-process supervisor: wait for the child, spool its output and serve followers.

The client owns all RunRecord state. This process owns only the child PID,
open pipes, transient FOLLOW socket and an acknowledged terminal receipt.
"""
from __future__ import annotations

import asyncio
import json
import os
import subprocess
import sys

from .remote_mailbox import publish_outcome
from .base import FailurePhase, RunState
from .remote_supervisor_stream import SupervisorStreams

_SPEC_ENV = "ASE2SPRKKR_REMOTE_PROCESS_SPEC"
_START_FD_ENV = "ASE2SPRKKR_REMOTE_PROCESS_START_FD"
_LAUNCH_FAILURE_RETURN_CODE = 127


def _write_outcome(receipt_id, state, return_code, failure_phase=None):
    publish_outcome(receipt_id, state, return_code, failure_phase)


def _wait_for_parent_release():
    value = os.environ.pop(_START_FD_ENV, None)
    if value is None:
        return
    fd = int(value)
    try:
        if os.read(fd, 1) != b"1":
            raise RuntimeError("Remote supervisor launch was not released by parent")
    finally:
        os.close(fd)


async def _pump(source, kind, streams):
    while True:
        data = await source.read(65536)
        if not data:
            return
        streams.publish(kind, data)


async def _run(spec):
    streams = SupervisorStreams(spec["output_file"], spec["error_file"])
    socket_path = await streams.start()
    # The parent reads this once before returning the backend handle; never
    # advertise a socket that is not listening yet.
    os.write(int(spec["socket_fd"]), (socket_path + "\n").encode())
    os.close(int(spec["socket_fd"]))
    opened = []
    stdin = subprocess.DEVNULL
    return_code = _LAUNCH_FAILURE_RETURN_CODE
    state = RunState.FAILED
    phase = FailurePhase.LAUNCH
    try:
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)
        try:
            process = await asyncio.create_subprocess_exec(
                *spec["argv"], cwd=spec["directory"], stdin=stdin,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT if streams.merged else asyncio.subprocess.PIPE,
            )
        except OSError:
            pass
        else:
            pumps = [asyncio.create_task(_pump(process.stdout, "out", streams))]
            if not streams.merged:
                pumps.append(asyncio.create_task(_pump(process.stderr, "err", streams)))
            return_code = await process.wait()
            await asyncio.gather(*pumps)
            state = RunState.FINISHED if return_code == 0 else RunState.FAILED
            phase = None if return_code == 0 else FailurePhase.EXECUTION
    finally:
        for file in opened:
            file.close()
        try:
            _write_outcome(spec["receipt_id"], state, return_code, phase)
            streams.terminal({
                "state": state.value, "return_code": return_code,
                "failure_phase": None if phase is None else phase.value,
            })
            await streams.finish()
        finally:
            pass
    return return_code


def main():
    raw = os.environ.pop(_SPEC_ENV, None)
    if raw is None:
        raise RuntimeError(f"Missing {_SPEC_ENV}")
    spec = json.loads(raw)
    _wait_for_parent_release()
    return asyncio.run(_run(spec))


if __name__ == "__main__":
    sys.exit(main())
