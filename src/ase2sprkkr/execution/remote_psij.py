"""PSI/J backend used by the detached remote execution helper."""

from __future__ import annotations

import asyncio
import os
import signal
import time
from dataclasses import asdict

from .base import ExecutionSpec, RunReference, RunState
from .psij import PsijExecutor


_FINAL_RUN_STATES = {
    RunState.FINISHED,
    RunState.FAILED,
    RunState.CANCELLED,
}
_FINAL_STATE_NAMES = {state.value for state in _FINAL_RUN_STATES}


class _RemoteResult:
    return_code = None


class _RemoteReader:
    """Drain scheduler output; parsing is done by the local RemoteRunner."""

    print_output = False

    def create_result(self, output_file, input_file):
        return _RemoteResult()

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.return_code = return_code
        return result


def _executor(request):
    options = request.get("executor_options") or {}
    return PsijExecutor(
        options.get("backend") or "slurm",
        url=options.get("url"),
    )


def run_remote_psij_worker(directory, run_id, request):
    """Run one remote request through PSI/J on the target machine."""
    from .remote import (
        _read_state,
        _read_state_unlocked,
        _state_lock,
        _timestamp,
        _write_state_unlocked,
    )

    output_file = request.get("output_file") or os.path.join(directory, "stdout")
    error_file = request.get("error_file") or output_file + ".err"
    spec = ExecutionSpec(
        argv=tuple(request["argv"]),
        input_file=request.get("input_file"),
        output_file=output_file,
        error_file=error_file,
        directory=os.fspath(directory),
        stdin=request.get("stdin"),
    )
    runner = _executor(request).create_run(_RemoteReader(), spec)

    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATE_NAMES:
            return
        if state.get("state") != RunState.QUEUED.value:
            raise RuntimeError(
                f"Cannot submit remote PSI/J run {run_id!r} from state "
                f"{state.get('state')!r}"
            )

        runner.submit()
        reference = runner.reference
        if reference is None:
            raise RuntimeError("PSI/J did not provide a persistent job reference")

        current = runner.state
        _write_state_unlocked(
            directory,
            run_id,
            state=current.value,
            pid=os.getpid(),
            backend_reference=asdict(reference),
            started_at=_timestamp(),
        )

    while current not in _FINAL_RUN_STATES:
        state = _read_state(directory, run_id)
        if state.get("state") == RunState.CANCELLED.value:
            runner.cancel()
        time.sleep(0.1)
        next_state = runner.state
        if next_state != current:
            current = next_state
            with _state_lock(directory, run_id):
                state = _read_state_unlocked(directory, run_id)
                if state.get("state") != RunState.CANCELLED.value:
                    _write_state_unlocked(directory, run_id, state=current.value)

    result = runner.wait()
    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        final_state = state.get("state")
        if final_state != RunState.CANCELLED.value:
            final_state = runner.state.value
        _write_state_unlocked(
            directory,
            run_id,
            state=final_state,
            return_code=result.return_code,
            finished_at=state.get("finished_at") or _timestamp(),
        )


def cancel_remote_psij(directory, run_id, request, timeout=1.0):
    """Cancel the remote PSI/J job, then stop its detached helper process."""
    from .remote import (
        _read_state,
        _read_state_unlocked,
        _state_lock,
        _timestamp,
        _write_state_unlocked,
    )

    with _state_lock(directory, run_id):
        state = _read_state_unlocked(directory, run_id)
        if state.get("state") in _FINAL_STATE_NAMES:
            return state

        reference = state.get("backend_reference")
        state = _write_state_unlocked(
            directory,
            run_id,
            state=RunState.CANCELLED.value,
            finished_at=_timestamp(),
        )
        pid = state.get("pid")

    if reference is not None:
        runner = _executor(request).attach(
            RunReference(**reference),
            _RemoteReader(),
        )
        runner.cancel()

    if not pid:
        return state

    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return _read_state(directory, run_id)

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return _read_state(directory, run_id)
        time.sleep(0.05)

    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return _read_state(directory, run_id)
