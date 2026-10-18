"""Direct-process supervisor: wait for the child, spool its output and serve followers.

The client owns all RunRecord state. This process owns only the child PID,
open pipes, transient FOLLOW socket and an acknowledged terminal receipt.
"""
from __future__ import annotations

import asyncio
import json
import os
import signal
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
    streams = SupervisorStreams(spec["output_file"], spec["error_file"], detached=spec.get("detached", True))
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
    process = None
    try:
        if not streams.detached:
            try:
                await asyncio.wait_for(streams.connected.wait(), timeout=30)
            except asyncio.TimeoutError:
                # No client took ownership of this non-persistent submission.
                return _LAUNCH_FAILURE_RETURN_CODE
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)
        try:
            process = await asyncio.create_subprocess_exec(
                *spec["argv"], cwd=spec["directory"], stdin=stdin,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT if streams.merged else asyncio.subprocess.PIPE,
                start_new_session=True,
            )
        except OSError:
            pass
        else:
            pumps = [asyncio.create_task(_pump(process.stdout, "out", streams))]
            if not streams.merged:
                pumps.append(asyncio.create_task(_pump(process.stderr, "err", streams)))
            if streams.detached:
                return_code = await process.wait()
            else:
                finished = asyncio.create_task(process.wait())
                lost = asyncio.create_task(streams.disconnected.wait())
                done, pending = await asyncio.wait((finished, lost), return_when=asyncio.FIRST_COMPLETED)
                if lost in done and not finished.done():
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        await asyncio.wait_for(asyncio.shield(finished), timeout=1)
                    except asyncio.TimeoutError:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        await finished
                    state = RunState.CANCELLED
                    phase = None
                for task in pending:
                    task.cancel()
                return_code = finished.result()
            await asyncio.gather(*pumps)
            if state != RunState.CANCELLED:
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
        finally:
            await streams.finish()
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
