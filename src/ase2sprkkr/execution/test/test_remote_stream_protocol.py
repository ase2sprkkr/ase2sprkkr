import asyncio
import sys
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec
from ..remote_endpoint import RemoteEndpointExecutor
from ..stream_protocol import (
    STREAM_MAGIC,
    FramedRunningCommand,
    StreamProtocolError,
    encode_stream_frame,
)
from ..transport import LocalTransport
from .test_remote import _launch, _processed_run


def _reader_with(data):
    stream = asyncio.StreamReader()
    stream.feed_data(data)
    stream.feed_eof()
    return stream


class _FinishedCommand:
    def __init__(self, stdout, stderr=b"", returncode=0):
        self.stdout = _reader_with(stdout)
        self.stderr = _reader_with(stderr)
        self.returncode = None
        self._returncode = returncode
        self.terminated = False

    async def wait(self):
        self.returncode = self._returncode
        return self.returncode

    def terminate(self):
        self.terminated = True


def test_framed_command_separates_calculation_and_transport_stderr():
    async def check():
        process = _FinishedCommand(
            STREAM_MAGIC
            + encode_stream_frame("out", b"stdout\n")
            + encode_stream_frame("err", b"stderr\n"),
            b"Warning: transport diagnostic\n",
        )
        command = FramedRunningCommand(process)

        stdout, stderr, returncode = await asyncio.gather(
            command.stdout.read(),
            command.stderr.read(),
            command.wait(),
        )

        assert returncode == 0
        assert stdout == b"stdout\n"
        assert stderr == b"stderr\n"
        assert command.diagnostics == b"Warning: transport diagnostic\n"

    asyncio.run(check())


def test_framed_command_reports_helper_failure_with_diagnostic():
    async def check():
        process = _FinishedCommand(
            STREAM_MAGIC,
            b"remote helper traceback\n",
            returncode=2,
        )
        command = FramedRunningCommand(process)

        with pytest.raises(StreamProtocolError, match="status 2.*remote helper traceback"):
            await command.wait()

    asyncio.run(check())


class _DiagnosticTransport(LocalTransport):
    """Inject an SSH-like diagnostic into the transport control channel."""

    async def start(self, argv, *, cwd=None, stdin=None):
        process = await super().start(argv, cwd=cwd, stdin=stdin)
        diagnostic = asyncio.StreamReader()
        diagnostic.feed_data(b"Warning: Permanently added test host\n")
        diagnostic.feed_eof()
        return SimpleNamespace(
            stdout=process.stdout,
            stderr=diagnostic,
            wait=process.wait,
            terminate=process.terminate,
            returncode=process.returncode,
        )


def test_remote_transport_stderr_never_enters_calculation_files(tmp_path):
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    input_file = (workdir / "input.txt").open("w+b")
    output_file = (workdir / "output.txt").open("wb")
    error_file = (workdir / "output.err").open("wb")
    remote_root = tmp_path / "remote"

    executor = RemoteEndpointExecutor(
        _DiagnosticTransport(),
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; print('stdout'); print('stderr', file=sys.stderr)",
        ),
        input_file=input_file,
        output_file=output_file,
        error_file=error_file,
        directory=str(workdir),
        launch=_launch(local_root, remote_root),
    )

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert result.output == b"stdout\n"
    assert result.error == b"stderr\n"
    assert (workdir / "output.txt").read_bytes() == b"stdout\n"
    assert (workdir / "output.err").read_bytes() == b"stderr\n"
    assert b"Permanently added" not in (workdir / "output.txt").read_bytes()
    assert b"Permanently added" not in (workdir / "output.err").read_bytes()


def test_framed_command_delivers_terminal_result_in_band():
    from ..stream_protocol import encode_terminal_frame

    async def check():
        command = FramedRunningCommand(_FinishedCommand(
            STREAM_MAGIC + encode_stream_frame("out", b"hello")
            + encode_terminal_frame({
                "state": "failed", "return_code": 4,
                "failure_phase": "execution", "error": "calculation failed",
            })
        ))
        assert await command.stdout.read() == b"hello"
        assert await command.wait() == 0
        assert command.terminal["return_code"] == 4
        assert command.terminal["failure_phase"] == "execution"
    asyncio.run(check())
