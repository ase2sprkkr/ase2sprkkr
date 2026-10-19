"""Output handle lifecycle shared by local, persistent and remote runners."""

import asyncio

import pytest

from ase2sprkkr.execution.base import ExecutionSpec, RunState, close_output_files
from ase2sprkkr.execution.local import LocalRunner
from ase2sprkkr.execution.local_project import LocalProjectRunner
from ase2sprkkr.execution.remote import RemoteRunner


class _TrackedOutput:
    def __init__(self, name=None):
        self.name = name
        self.close_calls = 0

    def close(self):
        self.close_calls += 1

    def flush(self):
        pass


@pytest.mark.parametrize("backend", ("local", "local_project"))
@pytest.mark.parametrize("merged", (False, True))
def test_local_runners_close_unique_output_handles(monkeypatch, tmp_path, backend, merged):
    stdout = _TrackedOutput(str(tmp_path / "stdout.out"))
    stderr = stdout if merged else _TrackedOutput(str(tmp_path / "stderr.err"))
    spec = ExecutionSpec(
        argv=("true",), directory=str(tmp_path),
        input_file=None, output_file=stdout, error_file=stderr,
    )

    if backend == "local":
        async def fail_to_launch(_spec):
            raise RuntimeError("launch failed")

        monkeypatch.setattr("ase2sprkkr.execution.local.start_local_process", fail_to_launch)
        runner = LocalRunner(spec)
        with pytest.raises(RuntimeError, match="launch failed"):
            asyncio.run(runner._execute())
        assert runner.state == RunState.FAILED
    else:
        runner = LocalProjectRunner(None, spec)
        monkeypatch.setattr(runner, "_monitor", lambda: None)
        runner._thread_main()
        assert runner._done.is_set()
        assert runner.exception is None

    assert stdout.close_calls == 1
    assert stderr.close_calls == 1


def test_close_output_files_skips_missing_and_path_declarations():
    first, second = _TrackedOutput(), _TrackedOutput()
    close_output_files(None, "out.txt", first, first, second)
    assert (first.close_calls, second.close_calls) == (1, 1)


def test_close_error_does_not_change_previous_failure_behavior():
    class Faulty:
        def close(self):
            raise OSError("cannot close output")

    other = _TrackedOutput()
    with pytest.raises(OSError, match="cannot close output"):
        close_output_files(Faulty(), other)
    assert other.close_calls == 0


def test_remote_monitor_closes_merged_output_once_after_follow_failure(tmp_path):
    output_path = tmp_path / "output.out"
    output_path.touch()
    output = _TrackedOutput(str(output_path))

    class UnreachableEndpoint:
        async def _follow(self, directory, backend_id, files, *, offsets):
            assert backend_id == "job-1"
            raise ConnectionError("remote follow interrupted")

    spec = ExecutionSpec(
        argv=("true",), directory=str(tmp_path),
        input_file=None, output_file=output, error_file=output,
    )
    runner = RemoteRunner(
        UnreachableEndpoint(), spec,
        backend_id="job-1", remote_directory="/remote/case",
    )
    asyncio.run(runner._monitor())

    assert runner.state == RunState.UNKNOWN
    assert output.close_calls == 1
