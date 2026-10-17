import asyncio

from ase2sprkkr.execution import ExecutionSpec, RunState, FileManifest, FileEntry
from ase2sprkkr.execution.remote import FileSyncPolicy, RemoteRunner
from ase2sprkkr.execution.transport import AtomicLocalFile


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

    async def _follow(self, directory, run_id, reference=None):
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
    output = AtomicLocalFile(destination)
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
        files=FileManifest(stdout=FileEntry(str(destination), "/work/case/calculation.out")),
        state=RunState.RUNNING,
    )
    return runner, executor, destination, output


def test_terminal_remote_state_may_publish_streamed_output(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.FINISHED)

    result = asyncio.run(runner._monitor())

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert destination.read_bytes() == b"streamed\n"
    assert output._finished is True
    assert executor.sync_calls == 2


def test_unknown_remote_state_does_not_publish_or_download_output(tmp_path):
    runner, executor, destination, output = _run(tmp_path, RunState.UNKNOWN)

    asyncio.run(runner._monitor())

    assert runner.state == RunState.UNKNOWN
    assert runner.failure_phase is None
    assert destination.read_bytes() == b"old\n"
    assert output._finished is True
    assert not output.temporary_path.exists()
    assert executor.sync_calls == 0
