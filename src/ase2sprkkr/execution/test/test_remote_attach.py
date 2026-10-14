import asyncio
import sys
import time
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, LaunchSpec, RunReference, RunState
from ..remote_endpoint import RemoteEndpointExecutor
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


def _executor(*, backend="local"):
    return RemoteEndpointExecutor(
        LocalTransport(),
        backend=backend,
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )


def _launch(local_root, remote_root):
    return LaunchSpec(local_root=local_root, remote_root=remote_root)


def _wait_for_state(runner, expected, timeout=5):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if runner.state == expected:
            return
        time.sleep(0.01)
    pytest.fail(f"remote run did not reach {expected.value}; current state={runner.state.value}")


def test_remote_executor_attach_reconnects_running_run_and_restores_stage_out(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    input_file = (workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (workdir / "run.out").open("wb")
    generated = workdir / "generated.dat"
    cache = workdir / "cache.dat"
    cache.write_text("before")

    first_executor = _executor()
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            (
                "import pathlib,time; "
                "release=pathlib.Path('release'); "
                "while not release.exists(): time.sleep(0.01); "
                "pathlib.Path('generated.dat').write_text('generated'); "
                "pathlib.Path('cache.dat').write_text('after'); "
                "print('done')"
            ),
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file, generated),
        cache_files=(cache,),
        launch=_launch(local_root, remote_root),
    )
    original = first_executor.create_run(_Reader(), spec)
    original.submit()
    _wait_for_state(original, RunState.RUNNING)

    reference = original.reference
    assert reference is not None

    # A fresh endpoint object must be able to reconstruct all per-run staging
    # policy from the persisted request; no live runner state is required.
    attached = _executor().attach(reference, _Reader())
    (remote_root / "calculation" / "release").write_text("go")

    assert attached.reference == reference
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"done\n"
    assert (workdir / "run.out").read_bytes() == b"done\n"
    assert (workdir / "generated.dat").read_text() == "generated"
    assert (workdir / "cache.dat").read_text() == "after"

    assert original.wait().return_code == 0


def test_remote_executor_attach_replays_finished_run(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    input_file = (workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (workdir / "run.out").open("wb")

    first_executor = _executor()
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", "print('finished')"),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file,),
        launch=_launch(local_root, remote_root),
    )
    original = first_executor.create_run(_Reader(), spec)
    assert original.run().return_code == 0
    reference = original.reference
    assert reference is not None

    attached = _executor().attach(reference, _Reader())
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"finished\n"
    assert (workdir / "run.out").read_bytes() == b"finished\n"


def test_remote_executor_attach_rejects_reference_for_other_backend(tmp_path):
    executor = _executor()
    reference = RunReference(
        executor="remote:psij",
        native_id="job",
        directory="/remote/calculation",
    )

    with pytest.raises(ValueError, match="belongs to"):
        executor.attach(reference, _Reader())
