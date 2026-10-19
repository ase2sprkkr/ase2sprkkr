"""Tests for command/file transports."""

import asyncio
import hashlib
import pytest
import os
import sys

from ..transport import LocalTransport, SshTransport
from ...common.file_versions import FileVersion, file_digest


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


def test_local_transport_preserves_file_signature(tmp_path):
    async def run_test():
        transport = LocalTransport()
        source = tmp_path / "source"
        destination = tmp_path / "destination"
        source.write_bytes(b"content")
        stat = source.stat()
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns - 1_000_000_000))

        source_signature = await transport.file_signature(source)
        await transport.get_file(source, destination)

        assert await transport.file_signature(destination) == source_signature

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

def test_transport_and_file_version_use_identical_digest(tmp_path):
    source = tmp_path / "large.bin"
    payload = b"one" * (400_000) + b"two"
    source.write_bytes(payload)
    expected = hashlib.sha256(payload).hexdigest()

    assert file_digest(source) == expected
    assert file_digest(source, missing_ok=True) == expected
    assert FileVersion.capture(source, require_digest=True).digest == expected
    assert asyncio.run(LocalTransport().file_digest(source)) == expected


def test_transport_missing_file_is_optional_but_provenance_digest_is_strict(tmp_path):
    missing = tmp_path / "missing"
    assert file_digest(missing, missing_ok=True) is None
    assert asyncio.run(LocalTransport().file_digest(missing)) is None
    with pytest.raises(FileNotFoundError):
        file_digest(missing)
