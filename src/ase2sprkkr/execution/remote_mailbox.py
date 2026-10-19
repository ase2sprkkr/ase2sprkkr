"""Durable mailbox for terminal remote-process outcomes, not RunRecord state."""

import json
import os
import re
from pathlib import Path
from uuid import uuid4

from ._durability import durable_mkdir, durable_write_json, fsync_directory
from .base import FailurePhase, RunState

_RECEIPT = re.compile(r"[0-9a-f]{32}\Z")
_TERMINAL = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _path(receipt_id):
    if not isinstance(receipt_id, str) or not _RECEIPT.fullmatch(receipt_id):
        raise ValueError("Invalid remote receipt ID")
    root = os.environ.get("ASE2SPRKKR_REMOTE_STATE_DIR")
    if root:
        base = Path(root).expanduser()
    else:
        base = Path(os.environ.get("XDG_STATE_HOME", "~/.local/state")).expanduser() / "ase2sprkkr"
    return base / "remote-outcomes" / (receipt_id + ".json")


def publish_outcome(receipt_id, state, return_code, failure_phase=None, *, error=None):
    """Publish once: cancellation must not overwrite natural completion."""
    state = RunState(state)
    if state not in _TERMINAL:
        raise ValueError("Backend receipt requires a terminal state")
    path = _path(receipt_id)
    durable_mkdir(path.parent, parents=True, exist_ok=True, mode=0o700)
    temporary = path.parent / ("." + receipt_id + "." + uuid4().hex + ".prepared")
    try:
        durable_write_json(
            temporary,
            {
                "state": state.value,
                "return_code": None if return_code is None else int(return_code),
                "failure_phase": None if failure_phase is None else FailurePhase(failure_phase).value,
                "error": error,
            },
            separators=(",", ":"),
            mode=0o600,
        )
        try:
            os.link(temporary, path)
        except FileExistsError:
            return False
        fsync_directory(path.parent)
        return True
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        fsync_directory(path.parent)


def read_outcome(receipt_id):
    """Read non-destructively so interrupted clients can retry."""
    try:
        result = json.loads(_path(receipt_id).read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    if not isinstance(result, dict) or result.get("state") not in {s.value for s in _TERMINAL}:
        raise ValueError("Invalid terminal receipt")
    return result


def acknowledge_outcome(receipt_id):
    """Delete only after the client has durably committed its own run state."""
    path = _path(receipt_id)
    try:
        path.unlink()
    except FileNotFoundError:
        return False
    fsync_directory(path.parent)
    return True
