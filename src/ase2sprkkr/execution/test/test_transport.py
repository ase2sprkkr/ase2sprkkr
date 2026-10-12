"""Tests for command/file transports."""

import subprocess
import sys

from ..transport import LocalTransport, SshTransport


def test_local_transport_exec_stream_and_files(tmp_path):
    transport = LocalTransport()

    result = transport.exec((sys.executable, "-c", "print('ok')"))
    assert result.stdout == b"ok\n"

    process = transport.exec_stream(
        (sys.executable, "-c", "import sys; print('stream'); print('err', file=sys.stderr)")
    )
    stdout, stderr = process.communicate()
    assert process.returncode == 0
    assert stdout == b"stream\n"
    assert stderr == b"err\n"

    source = tmp_path / "source"
    copied = tmp_path / "subdir" / "copied"
    returned = tmp_path / "returned"
    source.write_bytes(b"content")

    transport.put_file(source, copied)
    assert transport.read_file(copied) == b"content"

    transport.get_file(copied, returned)
    assert returned.read_bytes() == b"content"


def test_ssh_transport_quotes_remote_command(monkeypatch):
    seen = {}

    def run(command, **kwargs):
        seen["command"] = command
        seen["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, b"out", b"")

    monkeypatch.setattr(subprocess, "run", run)

    transport = SshTransport("cluster", options=("-p", "2222"))
    result = transport.exec(("printf", "%s", "a b"), cwd="/tmp/a b")

    assert result.stdout == b"out"
    assert seen["command"] == [
        "ssh",
        "-p",
        "2222",
        "cluster",
        "cd '/tmp/a b' && exec printf %s 'a b'",
    ]
