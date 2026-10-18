"""Direct-process supervisor: the sole owner of child waitpid and cancellation."""
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
    while data := await source.read(65536):
        streams.publish(kind, data)


async def _stop_process(process, waiting):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        await asyncio.wait_for(asyncio.shield(waiting), 1.0)
    except asyncio.TimeoutError:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        await waiting


async def _run(spec):
    streams = SupervisorStreams(
        spec.get("output_file"), spec.get("error_file"),
        detached=spec.get("detached", True),
    )
    await streams.start()
    os.write(int(spec["socket_fd"]), (streams.socket_path + "\n").encode())
    os.close(int(spec["socket_fd"]))
    state, phase, return_code, message = RunState.FAILED, FailurePhase.LAUNCH, 127, None
    opened = []
    process = None
    try:
        if not streams.detached:
            try:
                await asyncio.wait_for(streams.connected.wait(), timeout=30)
            except asyncio.TimeoutError:
                message = "Client never connected to attached supervisor"
                return return_code
        stdin = subprocess.DEVNULL
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)
        process = await asyncio.create_subprocess_exec(
            *spec["argv"], cwd=spec["directory"], stdin=stdin,
            stdout=asyncio.subprocess.PIPE,
            stderr=(asyncio.subprocess.STDOUT if streams.merged
                    else asyncio.subprocess.PIPE),
            start_new_session=True,
        )
        state, phase, return_code = RunState.RUNNING, None, None
        waiting = asyncio.create_task(process.wait())
        pumps = [asyncio.create_task(_pump(process.stdout, "out", streams))]
        if not streams.merged:
            pumps.append(asyncio.create_task(_pump(process.stderr, "err", streams)))
        # A failing spool must terminate the process rather than deadlock a pipe.
        stop = asyncio.create_task(streams.cancelled.wait())
        lost = (asyncio.create_task(streams.disconnected.wait())
                if not streams.detached else None)
        try:
            while True:
                watches = [waiting, stop, *(p for p in pumps if not p.done())]
                if lost is not None:
                    watches.append(lost)
                done, _ = await asyncio.wait(watches, return_when=asyncio.FIRST_COMPLETED)
                failed = next((t.exception() for t in pumps if t.done() and not t.cancelled()
                               and t.exception() is not None), None)
                if failed is not None:
                    message = f"Failed to write process output: {failed}"
                    state, phase = RunState.FAILED, FailurePhase.EXECUTION
                    await _stop_process(process, waiting)
                    break
                if waiting.done():
                    break
                if stop in done or (lost is not None and lost in done):
                    state, phase = RunState.CANCELLED, None
                    await _stop_process(process, waiting)
                    break
            return_code = await waiting
            await asyncio.gather(*pumps)
            if state == RunState.RUNNING:
                state = RunState.FINISHED if return_code == 0 else RunState.FAILED
                phase = None if return_code == 0 else FailurePhase.EXECUTION
        finally:
            stop.cancel()
            if lost is not None:
                lost.cancel()
    except OSError as exc:
        message = str(exc)
        state, phase = RunState.FAILED, (FailurePhase.LAUNCH if process is None
                                         else FailurePhase.EXECUTION)
    except Exception as exc:
        message = str(exc)
        state, phase = RunState.FAILED, FailurePhase.EXECUTION
        if process is not None and process.returncode is None:
            await _stop_process(process, asyncio.create_task(process.wait()))
        return_code = process.returncode if process is not None else 127
    finally:
        for f in opened:
            f.close()
        try:
            streams.sync()
        except OSError as exc:
            state, phase = RunState.FAILED, FailurePhase.EXECUTION
            message = f"Failed to sync process output: {exc}"
        outcome = {"state": state.value, "return_code": return_code,
                   "failure_phase": None if phase is None else phase.value,
                   "error": message}
        try:
            if message is not None:
                publish_outcome(spec["receipt_id"], state, return_code, phase, error=message)
            else:
                _write_outcome(spec["receipt_id"], state, return_code, phase)
        finally:
            streams.terminal(outcome)
            await streams.finish()
    return return_code if return_code is not None else 127


def main():
    raw = os.environ.pop(_SPEC_ENV, None)
    if raw is None:
        raise RuntimeError(f"Missing {_SPEC_ENV}")
    spec = json.loads(raw)
    _wait_for_parent_release()
    return asyncio.run(_run(spec))


if __name__ == "__main__":
    sys.exit(main())
