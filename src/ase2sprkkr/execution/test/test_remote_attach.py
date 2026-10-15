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
                "import pathlib,time\n"
                "release=pathlib.Path('release')\n"
                "while not release.exists():\n"
                "    time.sleep(0.01)\n"
                "pathlib.Path('generated.dat').write_text('generated')\n"
                "pathlib.Path('cache.dat').write_text('after')\n"
                "print('done')\n"
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


def test_remote_output_integrity_accepts_merged_local_stderr(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    output_file = (workdir / "run.out").open("wb")
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; print('stdout'); print('stderr', file=sys.stderr)",
        ),
        input_file=None,
        output_file=output_file,
        error_file=None,
        directory=str(workdir),
        launch=_launch(local_root, remote_root),
    )

    result = _executor().create_run(_Reader(), spec).run()

    assert result.return_code == 0
    assert result.output == b"stdout\n"
    assert result.error == b"stderr\n"
    local_log = (workdir / "run.out").read_bytes()
    assert len(local_log) == len(result.output) + len(result.error)
    assert b"stdout\n" in local_log
    assert b"stderr\n" in local_log
    assert (remote_root / "calculation" / "run.out").read_bytes() == b"stdout\n"
    assert (remote_root / "calculation" / "run.out.err").read_bytes() == b"stderr\n"


def test_remote_output_integrity_with_separate_error_file(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    output_file = (workdir / "run.out").open("wb")
    error_file = (workdir / "run.err").open("wb")
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; print('stdout'); print('stderr', file=sys.stderr)",
        ),
        input_file=None,
        output_file=output_file,
        error_file=error_file,
        directory=str(workdir),
        launch=_launch(local_root, remote_root),
    )

    result = _executor().create_run(_Reader(), spec).run()

    assert result.return_code == 0
    assert result.output == b"stdout\n"
    assert result.error == b"stderr\n"
    assert (workdir / "run.out").read_bytes() == b"stdout\n"
    assert (workdir / "run.err").read_bytes() == b"stderr\n"
    assert (remote_root / "calculation" / "run.out").read_bytes() == b"stdout\n"
    assert (remote_root / "calculation" / "run.err").read_bytes() == b"stderr\n"


def test_remote_executor_attach_rejects_reference_for_other_backend(tmp_path):
    executor = _executor()
    reference = RunReference(
        executor="remote:psij",
        native_id="job",
        directory="/remote/calculation",
    )

    with pytest.raises(ValueError, match="belongs to"):
        executor.attach(reference, _Reader())
