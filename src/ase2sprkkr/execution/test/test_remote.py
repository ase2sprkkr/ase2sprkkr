import asyncio
import sys
from types import SimpleNamespace

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


def test_remote_executor_local_transport_streams_output(tmp_path):
    input_file = (tmp_path / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    output_file = (tmp_path / "output.txt").open("wb")

    executor = RemoteExecutor(
        LocalTransport(),
        root_directory=str(tmp_path / "remote"),
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )
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
