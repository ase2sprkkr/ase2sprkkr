import asyncio
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, RunState
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


def _executor(local_root, remote_root):
    return RemoteExecutor(
        LocalTransport(),
        local_root=local_root,
        remote_root=remote_root,
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )


def test_remote_executor_local_transport_streams_output(tmp_path):
    input_file = (tmp_path / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    output_file = (tmp_path / "output.txt").open("wb")

    remote_root = tmp_path / "remote"
    executor = _executor(tmp_path, remote_root)
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; print(sys.stdin.buffer.read().decode().upper())",
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
    )

    seen = []
    runner = executor.create_run(_Reader(), spec)
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"HELLO\n"
    assert result.error == b""
    assert seen == [("HELLO\n", "out")]
    assert runner.reference is not None
    assert runner.reference.executor == "remote:local"
    assert Path(runner.reference.directory) == remote_root


def test_remote_directory_mapping_defaults_to_home():
    executor = RemoteExecutor(LocalTransport())
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == (Path.home() / "calculations" / "Fe").as_posix()


def test_remote_directory_mapping_can_replace_remote_root(tmp_path):
    executor = RemoteExecutor(
        LocalTransport(), remote_root="/scratch/tester"
    )
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == "/scratch/tester/calculations/Fe"


def test_remote_directory_mapping_can_replace_both_roots(tmp_path):
    local_root = tmp_path / "local"
    local = local_root / "project" / "Fe"
    executor = RemoteExecutor(
        LocalTransport(),
        local_root=local_root,
        remote_root="/scratch/tester/projects",
    )

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == "/scratch/tester/projects/project/Fe"


def test_remote_directory_mapping_rejects_path_outside_local_root(tmp_path):
    local_root = tmp_path / "local"
    executor = RemoteExecutor(
        LocalTransport(),
        local_root=local_root,
        remote_root="/scratch/tester",
    )

    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(executor.map_directory(tmp_path / "elsewhere"))
