import asyncio
import sys
from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionSpec, LocalExecutor, RunState
from ase2sprkkr.execution.base import LaunchSpec
from ase2sprkkr.execution.psij import PsijExecutor
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor
from ase2sprkkr.execution.transport import LocalTransport


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(stdout=None, stderr=None, return_code=None)

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.stdout = output
        result.stderr = error
        result.return_code = return_code
        return result


def _command():
    return (
        sys.executable,
        "-c",
        "import sys; "
        "sys.stdout.write('out\\n'); sys.stdout.flush(); "
        "sys.stderr.write('err\\n'); sys.stderr.flush()",
    )


def test_local_default_is_one_process_level_output_stream(tmp_path):
    output = (tmp_path / "run.out").open("wb")
    runner = LocalExecutor().create_run(
        _Reader(),
        ExecutionSpec(
            argv=_command(),
            input_file=None,
            output_file=output,
            directory=str(tmp_path),
        ),
    )
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))

    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert runner.reference.error_file is None
    assert result.return_code == 0
    assert result.stdout == b"out\nerr\n"
    assert result.stderr == b""
    assert (tmp_path / "run.out").read_bytes() == result.stdout
    assert {kind for kind, _ in seen} == {"out"}


def test_local_aliased_error_file_is_one_process_level_stream(tmp_path):
    output = (tmp_path / "run.out").open("wb")
    runner = LocalExecutor().create_run(
        _Reader(),
        ExecutionSpec(
            argv=_command(),
            input_file=None,
            output_file=output,
            error_file=output,
            directory=str(tmp_path),
        ),
    )

    result = runner.run()

    assert result.stdout == b"out\nerr\n"
    assert result.stderr == b""
    assert (tmp_path / "run.out").read_bytes() == result.stdout


def test_local_explicit_error_file_keeps_streams_separate(tmp_path):
    output = (tmp_path / "run.out").open("wb")
    error = (tmp_path / "run.err").open("wb")
    runner = LocalExecutor().create_run(
        _Reader(),
        ExecutionSpec(
            argv=_command(),
            input_file=None,
            output_file=output,
            error_file=error,
            directory=str(tmp_path),
        ),
    )

    result = runner.run()

    assert result.stdout == b"out\n"
    assert result.stderr == b"err\n"
    assert (tmp_path / "run.out").read_bytes() == b"out\n"
    assert (tmp_path / "run.err").read_bytes() == b"err\n"


def test_psij_recovery_hides_internal_merged_stderr_sidecar(tmp_path):
    output = str(tmp_path / "run.out")
    reference = PsijExecutor._recovered_reference(
        "slurm",
        None,
        "123",
        str(tmp_path),
        output,
        output + ".ase2sprkkr-psij-stderr",
    )

    assert reference.output_file == output
    assert reference.error_file is None


def test_remote_default_uses_same_authoritative_output(tmp_path):
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    input_file = (workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (workdir / "run.out").open("wb")
    spec = ExecutionSpec(
        argv=_command(),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file,),
        launch=LaunchSpec(
            local_root=str(local_root),
            remote_root=str(tmp_path / "remote"),
        ),
    )
    executor = RemoteEndpointExecutor(
        LocalTransport(),
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )

    runner = executor.create_run(_Reader(), spec)
    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert runner.reference.error_file is None
    assert result.stdout == b"out\nerr\n"
    assert result.stderr == b""
    assert (workdir / "run.out").read_bytes() == result.stdout
