"""Minimal direct-process supervisor used by the stateless remote service.

It owns no ASE2SPRKKR lifecycle state. Its only durable metadata is the backend
outcome needed after the supervised process has disappeared.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

from ._durability import durable_write_json
from .base import FailurePhase, RunState


_SPEC_ENV = "ASE2SPRKKR_REMOTE_PROCESS_SPEC"
_START_FD_ENV = "ASE2SPRKKR_REMOTE_PROCESS_START_FD"
_LAUNCH_FAILURE_RETURN_CODE = 127


def _write_outcome(path, state, return_code, failure_phase=None):
    """Atomically and durably publish one immutable backend terminal outcome."""
    payload = {
        "state": RunState(state).value,
        "return_code": None if return_code is None else int(return_code),
        "failure_phase": None if failure_phase is None else FailurePhase(failure_phase).value,
    }
    durable_write_json(path, payload, separators=(",", ":"))


def _wait_for_parent_release():
    value = os.environ.pop(_START_FD_ENV, None)
    if value is None:
        return
    fd = int(value)
    try:
        os.read(fd, 1)
    finally:
        os.close(fd)


def main():
    raw = os.environ.pop(_SPEC_ENV, None)
    if raw is None:
        raise RuntimeError(f"Missing {_SPEC_ENV}")
    spec = json.loads(raw)

    _wait_for_parent_release()

    opened = []
    stdin = subprocess.DEVNULL
    return_code = _LAUNCH_FAILURE_RETURN_CODE
    state = RunState.FAILED
    failure_phase = FailurePhase.LAUNCH
    try:
        if spec.get("stdin"):
            stdin = open(spec["stdin"], "rb")
            opened.append(stdin)

        output = open(spec["output_file"], "wb")
        opened.append(output)
        if spec["error_file"] == spec["output_file"]:
            error = subprocess.STDOUT
        else:
            error = open(spec["error_file"], "wb")
            opened.append(error)

        try:
            process = subprocess.Popen(
                spec["argv"],
                cwd=spec["directory"],
                stdin=stdin,
                stdout=output,
                stderr=error,
                close_fds=True,
            )
        except OSError:
            return_code = _LAUNCH_FAILURE_RETURN_CODE
            state = RunState.FAILED
            failure_phase = FailurePhase.LAUNCH
        else:
            return_code = process.wait()
            state = RunState.FINISHED if return_code == 0 else RunState.FAILED
            failure_phase = None if return_code == 0 else FailurePhase.EXECUTION
    finally:
        for file in opened:
            file.close()
        _write_outcome(spec["outcome_file"], state, return_code, failure_phase)

    return return_code


if __name__ == "__main__":
    sys.exit(main())
