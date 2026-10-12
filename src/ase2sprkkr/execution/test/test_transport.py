"""Tests for command/file transports."""

import asyncio
import sys

from ..transport import LocalTransport, SshTransport


def test_local_transport_run_start_and_files(tmp_path):
    async def run_test():
        transport = LocalTransport()

        result = await transport.run((sys.executable, "-c", "print('ok')"))
        assert result.stdout == b"ok\n"

        command = await transport.start(
            (
                sys.executable,
                "-c",
                "import sys; print('stream'); print('err', file=sys.stderr)",
            )
        )
        stdout, stderr, returncode = await asyncio.gather(
            command.stdout.read(),
            command.stderr.read(),
            command.wait(),
        )
        assert returncode == 0
        assert stdout == b"stream\n"
        assert stderr == b"err\n"

        source = tmp_path / "source"
        copied = tmp_path / "subdir" / "copied"
        returned = tmp_path / "returned"
        source.write_bytes(b"content")

        await transport.put_file(source, copied)
        assert await transport.read_file(copied) == b"content"

        await transport.get_file(copied, returned)
        assert returned.read_bytes() == b"content"

    asyncio.run(run_test())


def test_local_transport_drains_stdout_and_stderr_concurrently():
    """Both pipes must be drained while the child is still running."""

    async def run_test():
        transport = LocalTransport()
        chunks = 64
        chunk_size = 8192
        command = await transport.start(
            (
                sys.executable,
                "-c",
                (
                    "import sys\n"
                    f"chunk = b'x' * {chunk_size}\n"
                    f"for _ in range({chunks}):\n"
                    "    sys.stdout.buffer.write(chunk)\n"
                    "    sys.stdout.buffer.flush()\n"
                    "    sys.stderr.buffer.write(chunk)\n"
                    "    sys.stderr.buffer.flush()\n"
                ),
            )
        )

        stdout, stderr, returncode = await asyncio.wait_for(
            asyncio.gather(
                command.stdout.read(),
                command.stderr.read(),
                command.wait(),
            ),
            timeout=10,
        )

        assert returncode == 0
        expected_size = chunks * chunk_size
        assert len(stdout) == expected_size
        assert len(stderr) == expected_size

    asyncio.run(run_test())


def test_ssh_transport_quotes_remote_command():
    transport = SshTransport("cluster", options=("-p", "2222"))

    assert transport._command(("printf", "%s", "a b"), cwd="/tmp/a b") == [
        "ssh",
        "-p",
        "2222",
        "cluster",
        "cd '/tmp/a b' && exec printf %s 'a b'",
    ]
