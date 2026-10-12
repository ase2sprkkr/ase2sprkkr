import asyncio
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, RunReference, RunState
from ..remote import RemoteExecutor
from ..transport import LocalTransport


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None, output=None, error=None)

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.output = output
        result.error = error
        result.return_code = return_code
        return result


def _executor(local_root, remote_root, *, backend="local"):
    return RemoteExecutor(
        LocalTransport(),
        backend=backend,
        local_root=local_root,
        remote_root=remote_root,
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )


def _wait_for_state(runner, expected, timeout=5):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if runner.state == expected:
            return
        time.sleep(0.01)
    pytest.fail(f"remote run did not reach {expected.value}; current state={runner.state.value}")


def test_remote_executor_attach_reconnects_running_run_and_restores_stage_out(tmp_path):
    first_root = tmp_path / "first"
    second_root = tmp_path / "second"
    remote_root = tmp_path / "remote"
    first_workdir = first_root / "calculation"
    first_workdir.mkdir(parents=True)

    input_file = (first_workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (first_workdir / "run.out").open("wb")
    generated = first_workdir / "generated.dat"
    cache = first_workdir / "cache.dat"
    cache.write_text("before")

    first_executor = _executor(first_root, remote_root)
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            (
                "import pathlib,time; "
                "time.sleep(0.25); "
                "pathlib.Path('generated.dat').write_text('generated'); "
                "pathlib.Path('cache.dat').write_text('after'); "
                "print('done')"
            ),
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(first_workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file, generated),
        cache_files=(cache,),
    )
    original = first_executor.create_run(_Reader(), spec)
    original.submit()
    _wait_for_state(original, RunState.RUNNING)

    reference = original.reference
    assert reference is not None

    second_executor = _executor(second_root, remote_root)
    attached = second_executor.attach(reference, _Reader())

    assert attached.reference == reference
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"done\n"
    assert (second_root / "calculation" / "run.out").read_bytes() == b"done\n"
    assert (second_root / "calculation" / "generated.dat").read_text() == "generated"
    assert (second_root / "calculation" / "cache.dat").read_text() == "after"

    # The original client may still be alive in this process, but reconnecting
    # must not depend on any of its in-memory state.
    assert original.wait().return_code == 0


def test_remote_executor_attach_replays_finished_run(tmp_path):
    first_root = tmp_path / "first"
    second_root = tmp_path / "second"
    remote_root = tmp_path / "remote"
    first_workdir = first_root / "calculation"
    first_workdir.mkdir(parents=True)

    input_file = (first_workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (first_workdir / "run.out").open("wb")

    first_executor = _executor(first_root, remote_root)
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", "print('finished')"),
        input_file=input_file,
        output_file=output_file,
        directory=str(first_workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file,),
    )
    original = first_executor.create_run(_Reader(), spec)
    assert original.run().return_code == 0
    reference = original.reference
    assert reference is not None

    attached = _executor(second_root, remote_root).attach(reference, _Reader())
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"finished\n"
    assert (second_root / "calculation" / "run.out").read_bytes() == b"finished\n"


def test_remote_executor_attach_rejects_reference_for_other_backend(tmp_path):
    executor = _executor(tmp_path / "local", tmp_path / "remote")
    reference = RunReference(
        executor="remote:psij",
        native_id="job",
        directory="/remote/calculation",
    )

    with pytest.raises(ValueError, match="belongs to"):
        executor.attach(reference, _Reader())
