"""Durability boundaries shared by detached local and PSI/J execution."""

import os

import pytest

from .. import detached
from ..base import FailurePhase


def _running_state(output_file, error_file):
    return {
        "id": "run",
        "state": "running",
        "output_file": os.fspath(output_file),
        "error_file": os.fspath(error_file),
        "failure_phase": None,
    }


def test_result_files_are_synced_before_detached_terminal_state(tmp_path, monkeypatch):
    output = tmp_path / "run.out"
    error = tmp_path / "run.err"
    output.write_text("output")
    error.write_text("error")
    monkeypatch.setattr(
        detached,
        "_read_state_unlocked",
        lambda workspace, run_id: _running_state(output, error),
    )

    events = []
    monkeypatch.setattr(
        detached,
        "sync_file",
        lambda path: events.append(("sync", os.fspath(path))),
    )
    monkeypatch.setattr(
        detached,
        "durable_write_json",
        lambda path, state: events.append(("state", state["state"])),
    )

    detached._write_state_unlocked(tmp_path, "run", state="finished")

    assert events == [
        ("sync", os.fspath(output)),
        ("sync", os.fspath(error)),
        ("state", "finished"),
    ]


def test_detached_launch_failure_does_not_commit_stale_result_files(
    tmp_path, monkeypatch
):
    output = tmp_path / "run.out"
    output.write_text("old output")
    monkeypatch.setattr(
        detached,
        "_read_state_unlocked",
        lambda workspace, run_id: _running_state(output, output),
    )

    events = []
    monkeypatch.setattr(
        detached,
        "sync_file",
        lambda path: events.append(("sync", os.fspath(path))),
    )
    monkeypatch.setattr(
        detached,
        "durable_write_json",
        lambda path, state: events.append(("state", state["state"])),
    )

    detached._write_state_unlocked(
        tmp_path,
        "run",
        state="failed",
        failure_phase=FailurePhase.LAUNCH.value,
    )

    assert events == [("state", "failed")]


def test_detached_finalization_failure_does_not_retry_result_sync(
    tmp_path, monkeypatch
):
    output = tmp_path / "run.out"
    output.write_text("output")
    monkeypatch.setattr(
        detached,
        "_read_state_unlocked",
        lambda workspace, run_id: _running_state(output, output),
    )

    events = []
    monkeypatch.setattr(
        detached,
        "sync_file",
        lambda path: events.append(("sync", os.fspath(path))),
    )
    monkeypatch.setattr(
        detached,
        "durable_write_json",
        lambda path, state: events.append(("state", state["failure_phase"])),
    )

    detached._write_state_unlocked(
        tmp_path,
        "run",
        state="failed",
        failure_phase=FailurePhase.FINALIZATION.value,
    )

    assert events == [("state", FailurePhase.FINALIZATION.value)]


def test_detached_result_sync_failure_records_finalization_failure(
    tmp_path, monkeypatch
):
    output = tmp_path / "run.out"
    output.write_text("output")
    monkeypatch.setattr(
        detached,
        "_read_state_unlocked",
        lambda workspace, run_id: _running_state(output, output),
    )

    def fail_sync(path):
        raise OSError("cannot synchronize result")

    writes = []
    monkeypatch.setattr(detached, "sync_file", fail_sync)
    monkeypatch.setattr(
        detached,
        "durable_write_json",
        lambda path, state: writes.append(dict(state)),
    )

    with pytest.raises(OSError, match="cannot synchronize result"):
        detached._write_state_unlocked(tmp_path, "run", state="finished")

    assert len(writes) == 1
    failure = writes[0]
    assert failure["state"] == "failed"
    assert failure["failure_phase"] == FailurePhase.FINALIZATION.value
    assert "cannot synchronize result" in failure["error"]
    assert failure["finished_at"] is not None
