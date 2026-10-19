import asyncio

import pytest

from ase2sprkkr.execution import ExecutionSpec, RunState, FileManifest, FileEntry
from ase2sprkkr.execution.remote import FileSyncPolicy, RemoteRunner


class _Follow:
    def __init__(self, output=b"streamed\n", error=b""):
        self.terminal = None
        self.stdout = asyncio.StreamReader()
        self.stderr = asyncio.StreamReader()
        self.stdout.feed_data(output)
        self.stdout.feed_eof()
        self.stderr.feed_data(error)
        self.stderr.feed_eof()

    async def wait(self):
        return 0


class _Executor:
    output_files = FileSyncPolicy.NEVER

    def __init__(self, state):
        self.backend_state = RunState(state)
        self.sync_calls = 0

    async def _follow(self, directory, run_id, reference=None, *, offsets=None):
        return _Follow()

    async def _status(self, directory, run_id, reference=None):
        return {
            "state": self.backend_state.value,
            "return_code": 0 if self.backend_state == RunState.FINISHED else None,
            "failure_phase": None,
        }


class _GateRunner(RemoteRunner):
    async def _sync_output_files(self, required):
        self.executor.sync_calls += 1
        await super()._sync_output_files(required)

    async def _sync_cache_files(self):
        self.executor.sync_calls += 1
        await super()._sync_cache_files()


def _run(tmp_path, state, *, original=b"old\n"):
    destination = tmp_path / "calculation.out"
    destination.write_bytes(original)
    output = destination.open("wb")
    executor = _Executor(state)
    runner = _GateRunner(
        executor,
        ExecutionSpec(
            argv=(),
            input_file=None,
            output_file=output,
            directory=str(tmp_path),
        ),
        backend_id="opaque-handle",
        remote_directory="/work/case",
        files=FileManifest(
            stdout=FileEntry(str(destination), "/work/case/calculation.out"),
        ).initialize_availability(),
        state=RunState.RUNNING,
    )
    return runner, executor, destination, output


def test_terminal_remote_state_may_publish_streamed_output(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.FINISHED)

    result = asyncio.run(runner._monitor())

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert destination.read_bytes() == b"streamed\n"
    assert output.closed is True
    assert executor.sync_calls == 2


def test_unknown_remote_state_does_not_publish_or_download_output(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.UNKNOWN)

    asyncio.run(runner._monitor())

    assert runner.state == RunState.UNKNOWN
    assert runner.failure_phase is None
    assert destination.read_bytes() == b"streamed\n"
    assert output.closed is True
    assert not list(tmp_path.glob(".*.tmp"))
    assert executor.sync_calls == 0


async def _unavailable_output(required):
    raise OSError("remote output could not be downloaded")


def test_successful_backend_with_failed_download_fails_run(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.FINISHED)
    runner._sync_output_files = _unavailable_output

    with pytest.raises(OSError, match="could not be downloaded"):
        asyncio.run(runner._monitor())

    assert runner.state == RunState.FAILED
    assert runner.failure_phase.value == "download"
    assert runner.return_code == 0
    assert runner.files.availability_state == "failed"
    assert "could not be downloaded" in runner.files.availability_message
    assert destination.read_bytes() == b"streamed\n"
    assert output.closed is True


def test_failed_backend_keeps_execution_failure_when_diagnostics_cannot_download(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.FAILED)
    runner._sync_output_files = _unavailable_output

    async def failed_status(*_args):
        return {
            "state": "failed",
            "return_code": 17,
            "failure_phase": "execution",
            "error": "SPR-KKR exited with status 17",
        }

    executor._status = failed_status
    result = asyncio.run(runner._monitor())

    assert runner.state == RunState.FAILED
    assert runner.failure_phase.value == "execution"
    assert runner.failure_message == "SPR-KKR exited with status 17"
    assert result.return_code == 17
    assert runner.files.availability_state == "failed"
    assert "could not be downloaded" in runner.files.availability_message
    assert destination.read_bytes() == b"streamed\n"
    assert output.closed is True


def test_completed_remote_download_marks_files_available(tmp_path):
    runner, _, destination, _ = _run(tmp_path, RunState.FINISHED)
    asyncio.run(runner._monitor())
    assert runner.state == RunState.FINISHED
    assert runner.files.availability_state == "available"
    assert runner.files.availability_message is None
    assert destination.read_bytes() == b"streamed\n"


def test_cancelled_backend_survives_diagnostic_download_failure(tmp_path):
    runner, _, destination, _ = _run(tmp_path, RunState.CANCELLED)
    runner._sync_output_files = _unavailable_output
    asyncio.run(runner._monitor())
    assert runner.state == RunState.CANCELLED
    assert runner.failure_phase is None
    assert runner.files.availability_state == "failed"
    assert destination.read_bytes() == b"streamed\n"



def test_remote_stderr_incomplete_after_success_is_download_failure(tmp_path):
    """FINISHED must not hide a truncated stderr FOLLOW stream."""
    from ase2sprkkr.execution.base import FailurePhase
    from ase2sprkkr.execution.remote_endpoint import _RemoteRun

    stdout_path = tmp_path / "run.out"
    stderr_path = tmp_path / "run.err"
    output = stdout_path.open("wb")
    error = stderr_path.open("wb")

    class Transport:
        async def file_signature(self, remote):
            return (len(b"streamed\n") if remote.endswith(".out") else 8, 0)

        async def file_digest(self, remote):
            import hashlib
            return hashlib.sha256(b"streamed\n").hexdigest()

    executor = _Executor(RunState.FINISHED)
    executor.transport = Transport()
    runner = _RemoteRun(
        executor,
        ExecutionSpec(
            argv=(), input_file=None, output_file=output,
            error_file=error, directory=str(tmp_path),
        ),
        backend_id="opaque-handle", remote_directory="/work/case",
        files=FileManifest(
            stdout=FileEntry(str(stdout_path), "/work/case/run.out"),
            stderr=FileEntry(str(stderr_path), "/work/case/run.err"),
        ).initialize_availability(),
        state=RunState.RUNNING,
    )

    with pytest.raises(RuntimeError, match="stderr size does not match"):
        asyncio.run(runner._monitor())

    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.DOWNLOAD
    assert runner.files.stderr.availability_state == "failed"
    assert runner.files.stdout.availability_state == "available"
    assert stdout_path.read_bytes() == b"streamed\n"
    assert stderr_path.read_bytes() == b""
    assert output.closed and error.closed



def test_follow_write_is_idempotent_across_attached_clients(tmp_path):
    """Two clients consuming the same chunks must not duplicate the output."""
    from ase2sprkkr.execution.remote import RemoteRunner
    destination = tmp_path / "output"
    destination.write_bytes(b"abc")
    first = destination.open("ab")
    second = destination.open("ab")
    try:
        spec1 = ExecutionSpec(argv=(), input_file=None, output_file=first, directory=str(tmp_path))
        spec2 = ExecutionSpec(argv=(), input_file=None, output_file=second, directory=str(tmp_path))
        a = RemoteRunner(_Executor(RunState.RUNNING), spec1)
        b = RemoteRunner(_Executor(RunState.RUNNING), spec2)
        a._stream_positions = {"out": 3}
        b._stream_positions = {"out": 3}
        a._write_output(b"def", "out")
        b._write_output(b"def", "out")
        assert destination.read_bytes() == b"abcdef"
        assert a._stream_positions["out"] == b._stream_positions["out"] == 6
        b._stream_positions["out"] = 3
        with pytest.raises(RuntimeError, match="Conflicting FOLLOW bytes"):
            b._write_output(b"wrong", "out")
        assert destination.read_bytes() == b"abcdef"
    finally:
        first.close()
        second.close()


def test_remote_follow_transport_failure_is_recoverable_unknown(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.RUNNING)

    async def disconnected(*_args, **_kwargs):
        raise ConnectionError("SSH disconnected")

    executor._follow = disconnected
    result = asyncio.run(runner._monitor())
    assert result.return_code is None
    assert runner._state == RunState.UNKNOWN
    assert runner.failure_phase is None
    assert runner.failure_message == "SSH disconnected"
    assert output.closed


def test_remote_status_transport_failure_after_follow_is_unknown(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.RUNNING)

    async def unavailable(*_args, **_kwargs):
        raise OSError("SSH lost during status query")

    executor._status = unavailable
    result = asyncio.run(runner._monitor())
    assert result.return_code is None
    assert runner._state == RunState.UNKNOWN
    assert runner.failure_phase is None
    assert runner.failure_message == "SSH lost during status query"
    assert output.closed
