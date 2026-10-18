import asyncio

import pytest

from ase2sprkkr.execution import ExecutionSpec, RunState, FileManifest, FileEntry
from ase2sprkkr.execution.remote import FileSyncPolicy, RemoteRunner


class _Follow:
    def __init__(self, output=b"streamed\n", error=b""):
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


async def _unavailable_output(_required):
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
