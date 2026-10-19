"""Regression tests for authoritative stdout/stderr verification after FOLLOW."""

import asyncio
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, FailurePhase, RunState
from ..remote import RemoteRunner
from ..remote_endpoint import _RemoteRun
from ..run_model import FileEntry, FileManifest
from ..transport import LocalTransport


async def _completed(self):
    self._state = RunState.FINISHED
    return self.result


def _make_runner(tmp_path, monkeypatch, *, stdout=b"same", stderr=b"error"):
    monkeypatch.setattr(RemoteRunner, "_monitor", _completed)
    local_out = tmp_path / "local.out"
    local_err = tmp_path / "local.err"
    remote_out = tmp_path / "remote.out"
    remote_err = tmp_path / "remote.err"
    local_out.write_bytes(stdout)
    local_err.write_bytes(stderr)
    remote_out.write_bytes(b"same")
    remote_err.write_bytes(b"error")
    runner = _RemoteRun(
        SimpleNamespace(transport=LocalTransport()),
        ExecutionSpec(argv=(), input_file=None, output_file=str(local_out),
                      error_file=str(local_err), directory=str(tmp_path)),
        backend_id="backend", remote_directory=str(tmp_path),
        files=FileManifest(
            stdout=FileEntry(str(local_out), str(remote_out)),
            stderr=FileEntry(str(local_err), str(remote_err)),
        ).initialize_availability(),
        state=RunState.FINISHED,
    )
    return runner, local_out, local_err, remote_out, remote_err


@pytest.mark.parametrize("channel", ("stdout", "stderr"))
def test_equal_size_corruption_fails_even_if_length_matches(tmp_path, monkeypatch, channel):
    runner, local_out, local_err, _, _ = _make_runner(
        tmp_path, monkeypatch,
        stdout=b"bad!" if channel == "stdout" else b"same",
        stderr=b"wrong" if channel == "stderr" else b"error",
    )
    with pytest.raises(RuntimeError, match=f"{channel} contents differ"):
        asyncio.run(runner._monitor())

    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.DOWNLOAD
    assert getattr(runner.files, channel).availability_state == "failed"
    assert (local_out if channel == "stdout" else local_err).exists()


def test_missing_authoritative_stream_is_not_success(tmp_path, monkeypatch):
    runner, _, _, _, remote_err = _make_runner(tmp_path, monkeypatch)
    remote_err.unlink()
    with pytest.raises(FileNotFoundError, match="Remote stderr output is missing"):
        asyncio.run(runner._monitor())
    assert runner.state == RunState.FAILED
    assert runner.files.stderr.availability_state == "failed"


def test_successful_integrity_check_preserves_timestamps(tmp_path, monkeypatch):
    runner, local_out, local_err, remote_out, remote_err = _make_runner(
        tmp_path, monkeypatch,
    )
    asyncio.run(runner._monitor())
    assert runner.state == RunState.FINISHED
    assert local_out.stat().st_mtime_ns == remote_out.stat().st_mtime_ns
    assert local_err.stat().st_mtime_ns == remote_err.stat().st_mtime_ns
