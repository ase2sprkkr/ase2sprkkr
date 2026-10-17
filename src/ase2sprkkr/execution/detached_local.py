"""Local process backend for the detached execution protocol."""

from __future__ import annotations

import asyncio
import os
from pathlib import Path

from ._durability import durable_mkdir, durable_replace, sync_stream
from .base import FailurePhase
from .local_process import (
    materialize_local_spec,
    start_local_process,
    stop_local_process,
)


def _launch_stream_path(path, run_id, kind):
    return path.with_name(f".{path.name}.{run_id}.{kind}.launching")


def _output_path(directory, value, default):
    path = Path(value or default)
    if not path.is_absolute():
        path = directory / path
    durable_mkdir(path.parent, parents=True, exist_ok=True)
    return path


async def _run_detached_local(directory, run_id, request):
    from .detached import (
        _FINAL_STATES,
        _read_state_unlocked,
        _request_spec,
        _state_lock,
        _timestamp,
        _write_state_unlocked,
    )

    spec = materialize_local_spec(_request_spec(directory, request))
    output_path = _output_path(directory, spec.output_file, "stdout")
    error_path = _output_path(directory, spec.error_file, "stderr")
    merged = error_path == output_path

    output_launch = _launch_stream_path(output_path, run_id, "out")
    error_launch = output_launch if merged else _launch_stream_path(
        error_path, run_id, "err"
    )
    launch_paths = (output_launch,) if merged else (output_launch, error_launch)

    if spec.stdin:
        stdin_path = Path(spec.stdin)
        if not stdin_path.is_absolute():
            stdin_path = directory / stdin_path
        stdin = stdin_path.open("rb")
    else:
        stdin = asyncio.subprocess.DEVNULL

    process = None
    return_code = None
    preserve_launch_paths = False

    def record_artifact_failure(exc):
        nonlocal preserve_launch_paths
        preserve_launch_paths = True
        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            if state.get("state") != "cancelled":
                _write_state_unlocked(
                    directory,
                    run_id,
                    state="failed",
                    return_code=return_code,
                    error=repr(exc),
                    failure_phase=FailurePhase.ARTIFACTS.value,
                    finished_at=state.get("finished_at") or _timestamp(),
                )

    try:
        with output_launch.open("wb") as output:
            if merged:
                error = None
                error_context = None
            else:
                error_context = error_launch.open("wb")
                error = error_context

            try:
                with _state_lock(directory, run_id):
                    state = _read_state_unlocked(directory, run_id)
                    if state.get("state") in _FINAL_STATES:
                        return
                    if state.get("state") != "submitting":
                        raise RuntimeError(
                            f"Cannot start detached run {run_id!r} from state "
                            f"{state.get('state')!r}"
                        )

                    process = await start_local_process(
                        spec,
                        stdin=stdin,
                        stdout=output,
                        stderr=(asyncio.subprocess.STDOUT if merged else error),
                        start_new_session=True,
                    )
                    try:
                        _write_state_unlocked(
                            directory,
                            run_id,
                            state="running",
                            pid=process.pid,
                            failure_phase=None,
                            stream_output_file=os.fspath(output_launch),
                            stream_error_file=os.fspath(error_launch),
                            started_at=_timestamp(),
                        )
                    except BaseException:
                        await stop_local_process(process, process_group=True)
                        raise

                try:
                    return_code = await process.wait()
                except BaseException:
                    await stop_local_process(process, process_group=True)
                    raise

                try:
                    sync_stream(output)
                    if error_context is not None:
                        sync_stream(error_context)
                except BaseException as exc:
                    record_artifact_failure(exc)
                    raise
            finally:
                if error_context is not None:
                    error_context.close()

        try:
            # Publication is part of artifact availability, not result parsing.
            # FINISHED is written only after these durable replacements succeed.
            durable_replace(output_launch, output_path)
            if not merged:
                durable_replace(error_launch, error_path)
        except BaseException as exc:
            record_artifact_failure(exc)
            raise

        with _state_lock(directory, run_id):
            state = _read_state_unlocked(directory, run_id)
            final_state = state.get("state")
            if final_state != "cancelled":
                final_state = "finished" if return_code == 0 else "failed"
            _write_state_unlocked(
                directory,
                run_id,
                state=final_state,
                output_file=os.fspath(output_path),
                error_file=os.fspath(error_path),
                stream_output_file=None,
                stream_error_file=None,
                return_code=return_code,
                failure_phase=(
                    FailurePhase.EXECUTION.value if final_state == "failed" else None
                ),
                finished_at=state.get("finished_at") or _timestamp(),
            )
    finally:
        if stdin is not asyncio.subprocess.DEVNULL:
            stdin.close()
        if not preserve_launch_paths:
            for path in launch_paths:
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass


def run_detached_local(directory, run_id, request):
    """Run one local calculation owned by the detached helper."""
    directory = Path(directory).expanduser().resolve()
    asyncio.run(_run_detached_local(directory, run_id, request))
