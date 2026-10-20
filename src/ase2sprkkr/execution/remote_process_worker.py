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
from .local_process import pump_stream
from .remote_supervisor_stream import SupervisorStreams

_SPEC_ENV = "ASE2SPRKKR_REMOTE_PROCESS_SPEC"
_START_FD_ENV = "ASE2SPRKKR_REMOTE_PROCESS_START_FD"


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


async def _launch_managed_process(spec, stdin, stderr):
    """Start the watchdog before launching a child into its process group.

    The watchdog and calculation share the supervisor's session, not its
    process group. The ready handshake closes the parent-death race.
    """
    read_fd, write_fd = os.pipe()
    watchdog = None
    try:
        watchdog = await asyncio.create_subprocess_exec(
            sys.executable, "-m", "ase2sprkkr.execution.remote_process_watchdog",
            str(read_fd),
            stdin=subprocess.DEVNULL, stdout=asyncio.subprocess.PIPE,
            stderr=subprocess.DEVNULL, pass_fds=(read_fd,),

        )
        os.close(read_fd)
        read_fd = None
        ready = await asyncio.wait_for(watchdog.stdout.readexactly(1), 5.0)
        if ready != b"R":
            raise RuntimeError("Calculation watchdog failed to initialize")
        # A tiny exec trampoline joins the group using os.setpgid().
        # Unlike preexec_fn, this does not run Python inside a forked,
        # potentially multi-threaded supervisor.
        error_read, error_write = os.pipe()
        try:
            try:
                process = await asyncio.create_subprocess_exec(
                    sys.executable, "-m",
                    "ase2sprkkr.execution.remote_process_watchdog",
                    "--exec", str(watchdog.pid), str(error_write), *spec["argv"],
                    cwd=spec["directory"], stdin=stdin,
                    stdout=asyncio.subprocess.PIPE, stderr=stderr,
                    pass_fds=(error_write,),
                )
            finally:
                os.close(error_write)
            # The CLOEXEC descriptor closes when the requested program
            # starts; errors before exec are sent through the same pipe.
            error = await asyncio.to_thread(os.read, error_read, 4096)
            if error:
                await process.wait()
                errno_number, message = json.loads(error.decode("utf8"))
                raise OSError(errno_number, message)
        finally:
            os.close(error_read)
        return process, watchdog, write_fd
    except BaseException:
        os.close(write_fd)
        if watchdog is not None:
            await watchdog.wait()
        raise
    finally:
        if read_fd is not None:
            os.close(read_fd)


async def _run(spec):
    streams = SupervisorStreams(
        spec.get("output_file"), spec.get("error_file"),
        detached=spec.get("detached", True),
        socket_dir=spec.get("socket_dir"),
    )
    await streams.start()
    os.write(int(spec["socket_fd"]), (streams.socket_path + "\n").encode())
    os.close(int(spec["socket_fd"]))
    state, phase, return_code, message = RunState.FAILED, FailurePhase.LAUNCH, 127, None
    opened = []
    process = None
    watchdog = None
    watchdog_fd = None
    try:
        if not streams.detached:
            try:
                await asyncio.wait_for(streams.connected.wait(), timeout=30)
            except asyncio.TimeoutError:
                message = "Client never connected to attached supervisor"
                return return_code
        if streams.cancelled.is_set() or (not streams.detached and streams.disconnected.is_set()):
            state, phase, return_code = RunState.CANCELLED, None, -signal.SIGTERM
            return return_code
        stdin = subprocess.DEVNULL
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)
        process, watchdog, watchdog_fd = await _launch_managed_process(
            spec, stdin,
            asyncio.subprocess.STDOUT if streams.merged else asyncio.subprocess.PIPE,
        )
        state, phase, return_code = RunState.RUNNING, None, None
        waiting = asyncio.create_task(process.wait())
        pumps = [asyncio.create_task(
            pump_stream(process.stdout, lambda data: streams.publish("out", data))
        )]
        if not streams.merged:
            pumps.append(asyncio.create_task(
                pump_stream(process.stderr, lambda data: streams.publish("err", data))
            ))
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
                elif not waiting.done():
                    if stop in done or (lost is not None and lost in done):
                        state, phase = RunState.CANCELLED, None
                    else:
                        continue

                # Stop descendants holding output pipes, then drain remaining
                # output after the main process has terminated.
                os.close(watchdog_fd)
                watchdog_fd = None
                break
            return_code = await waiting
            # Inspect all pump outcomes without allowing their exceptions to
            # replace the diagnosed output-write failure and its context.
            pump_results = await asyncio.gather(*pumps, return_exceptions=True)
            pump_error = next(
                (value for value in pump_results if isinstance(value, BaseException)),
                None,
            )
            if pump_error is not None:
                state, phase = RunState.FAILED, FailurePhase.EXECUTION
                message = f"Failed to write process output: {pump_error}"
            elif state == RunState.RUNNING:
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
        # The watchdog will stop the entire group in the finally block.
        return_code = process.returncode if process is not None else 127
    finally:
        if watchdog_fd is not None:
            os.close(watchdog_fd)
        if watchdog is not None:
            await watchdog.wait()
        if process is not None and process.returncode is None:
            await process.wait()
        if process is not None and return_code is None:
            return_code = process.returncode
        for f in opened:
            f.close()
        try:
            streams.sync()
        except OSError as exc:
            state, phase = RunState.FAILED, FailurePhase.ARTIFACTS
            message = f"Failed to sync process output: {exc}"
        outcome = {"state": state.value, "return_code": return_code,
                   "failure_phase": None if phase is None else phase.value,
                   "error": message}
        try:
            if message is not None:
                publish_outcome(spec["receipt_id"], state, return_code, phase, error=message)
            else:
                _write_outcome(spec["receipt_id"], state, return_code, phase)
        except OSError as exc:
            # No durable terminal receipt exists. The client must preserve its
            # active run and retry recovery, never commit FINISHED from a socket.
            outcome = {
                "state": RunState.UNKNOWN.value,
                "return_code": None,
                "failure_phase": None,
                "error": f"Unable to persist terminal receipt: {exc}",
            }
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
