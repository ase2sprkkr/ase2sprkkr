"""Tests for command/file transports."""

import asyncio
import hashlib
import pytest
import os
import signal
import subprocess
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



def test_ssh_and_local_globstar_agree_on_hidden_entries(tmp_path):
    """Remote glob semantics match the local transport, including dotfiles."""
    from ..transport import _match_relative_path, expand_local_files

    for name in (
        "shown/one.dat",
        "shown/.private.dat",
        ".hidden/other.dat",
        "shown/.hidden/deep.dat",
        "shown/nested/two.dat",
        "shown/nested/.another.dat",
    ):
        destination = tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(name)

    names = tuple(
        str(path.relative_to(tmp_path))
        for path in tmp_path.rglob("*")
        if path.is_file()
    )
    for pattern in (
        "**/*.dat", "shown/**/*.dat", "**/.*.dat",
        "**/.hidden/*.dat", ".hidden/*.dat", "shown/one.dat",
    ):
        local = {
            str(os.path.relpath(path, tmp_path))
            for path in expand_local_files(pattern, tmp_path)
        }
        ssh = {name for name in names if _match_relative_path(name, pattern)}
        assert ssh == local, (pattern, ssh, local)


def test_strict_transfer_refuses_wildcard_symlink_outside_source_root(tmp_path):
    import pytest
    from ..transport import expand_local_candidates, expand_local_files

    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    target = outside / "sensitive.dat"
    target.write_text("protected")
    alias = root / "alias.dat"
    alias.symlink_to(target)

    assert expand_local_candidates(str(root / "*.dat")) == (str(alias),)
    with pytest.raises(ValueError, match="outside source root"):
        expand_local_files(str(root / "*.dat"), root)
    assert target.read_text() == "protected"



class _LoopbackSshTransport(SshTransport):
    """Exercise the real SSH transfer implementation without an SSH server."""

    def __init__(self, *, zstd=False):
        super().__init__("local-test", zstd=zstd)

    def _command(self, argv, cwd=None):
        assert cwd is None
        return tuple(os.fspath(part) for part in argv)


def test_ssh_batch_tar_upload_and_download_share_pipeline(tmp_path):
    import shutil

    if shutil.which("tar") is None:
        pytest.skip("tar executable unavailable")

    source = tmp_path / "upload source"
    (source / "nested").mkdir(parents=True)
    (source / "nested" / "one.dat").write_text("one")
    (source / "nested" / "two.dat").write_text("two")
    remote = tmp_path / "remote target"
    downloaded = tmp_path / "download destination"
    transport = _LoopbackSshTransport()

    async def check():
        sent = await transport.put_files("nested/*.dat", source, remote)
        assert len(sent) == 2
        assert all(os.path.isfile(destination) for _, destination in sent)
        received = await transport.get_files("nested/*.dat", remote, downloaded)
        assert len(received) == 2
        assert all(os.path.isfile(destination) for _, destination in received)

    asyncio.run(check())
    for name in ("one.dat", "two.dat"):
        expected = (source / "nested" / name).read_text()
        assert (remote / "nested" / name).read_text() == expected
        assert (downloaded / "nested" / name).read_text() == expected


@pytest.mark.parametrize("failing_side", ["producer", "consumer"])
def test_ssh_batch_tar_pipeline_reports_subprocess_failures(failing_side):
    transport = _LoopbackSshTransport()
    producer = (
        sys.executable, "-c",
        "import sys; sys.stdout.buffer.write(b'data'); "
        "sys.stdout.buffer.flush(); "
        "sys.stderr.write('producer failed'); "
        f"sys.exit({7 if failing_side == 'producer' else 0})",
    )
    consumer = (
        sys.executable, "-c",
        "import sys; sys.stdin.buffer.read(); "
        "sys.stderr.write('consumer failed'); "
        f"sys.exit({9 if failing_side == 'consumer' else 0})",
    )
    async def check():
        with pytest.raises(__import__("subprocess").CalledProcessError) as exc:
            await asyncio.wait_for(transport._transfer_tar(producer, consumer), timeout=5)
        expected = 7 if failing_side == "producer" else 9
        assert exc.value.returncode == expected
        assert (failing_side + " failed").encode() in exc.value.stderr

    asyncio.run(check())


def test_ssh_batch_tar_pipeline_drains_stderr_during_copy():
    transport = _LoopbackSshTransport()
    producer = (
        sys.executable, "-c",
        "import sys; sys.stderr.buffer.write(b'e' * 262144); "
        "sys.stderr.buffer.flush(); sys.stdout.buffer.write(b'payload')",
    )
    consumer = (
        sys.executable, "-c",
        "import sys; assert sys.stdin.buffer.read() == b'payload'",
    )
    asyncio.run(asyncio.wait_for(
        transport._transfer_tar(producer, consumer), timeout=5,
    ))


async def _wait_for_command_pid(path):
    for _ in range(500):
        if path.exists():
            return int(path.read_text())
        await asyncio.sleep(0.01)
    raise AssertionError("Subprocess never reported its PID")


async def _cancel_and_check_reaped(task, pid):
    task.cancel()
    with pytest.raises(asyncio.CancelledError):
        await asyncio.wait_for(task, timeout=5)
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


@pytest.mark.skipif(os.name != "posix", reason="PID cleanup assertion is POSIX-only")
@pytest.mark.parametrize("factory", (LocalTransport, _LoopbackSshTransport))
def test_transport_run_cancellation_reaps_subprocess(tmp_path, factory):
    pid_file = tmp_path / "command.pid"

    async def check():
        command = (
            sys.executable, "-c",
            f"import os, time; open({str(pid_file)!r}, 'w').write(str(os.getpid())); "
            "time.sleep(30)",
        )
        task = asyncio.create_task(factory().run(command))
        pid = None
        try:
            pid = await _wait_for_command_pid(pid_file)
            await _cancel_and_check_reaped(task, pid)
        finally:
            if not task.done():
                task.cancel()
                await asyncio.gather(task, return_exceptions=True)
            if pid is not None:
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    asyncio.run(check())


@pytest.mark.skipif(os.name != "posix", reason="Requires POSIX shell and PID semantics")
def test_ssh_get_file_cancellation_reaps_process_and_removes_temp(tmp_path):
    pid_file = tmp_path / "download.pid"
    source = tmp_path / "source"
    destination = tmp_path / "downloaded"
    source.write_bytes(b"source bytes")
    destination.write_bytes(b"old destination")

    class SlowDownloadTransport(_LoopbackSshTransport):
        def _command(self, argv, cwd=None):
            if argv[0] == "cat":
                return (
                    sys.executable, "-c",
                    f"import os, time; open({str(pid_file)!r}, 'w').write(str(os.getpid())); "
                    "time.sleep(30)",
                )
            return super()._command(argv, cwd)

    async def check():
        task = asyncio.create_task(SlowDownloadTransport().get_file(source, destination))
        pid = None
        try:
            pid = await _wait_for_command_pid(pid_file)
            await _cancel_and_check_reaped(task, pid)
        finally:
            if not task.done():
                task.cancel()
                await asyncio.gather(task, return_exceptions=True)
            if pid is not None:
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    asyncio.run(check())
    assert destination.read_bytes() == b"old destination"
    assert not tuple(tmp_path.glob(".downloaded.*.tmp"))


@pytest.mark.skipif(os.name != "posix", reason="Requires POSIX shell tools")
def test_ssh_get_file_preserves_subprocess_failure(tmp_path):
    source = tmp_path / "source"
    source.write_bytes(b"source bytes")
    destination = tmp_path / "downloaded"

    class FailingDownloadTransport(_LoopbackSshTransport):
        def _command(self, argv, cwd=None):
            if argv[0] == "cat":
                return (
                    sys.executable, "-c",
                    "import sys; sys.stderr.write('download failed'); sys.exit(17)",
                )
            return super()._command(argv, cwd)

    with pytest.raises(subprocess.CalledProcessError) as failure:
        asyncio.run(FailingDownloadTransport().get_file(source, destination))
    assert failure.value.returncode == 17
    assert b"download failed" in failure.value.stderr
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".downloaded.*.tmp"))


@pytest.mark.parametrize(
    "method, command, stdout, expected",
    [
        ("file_digest", 'sha256sum -- "$1"',
         b"1234abcd  folder name/result.out\n", "1234abcd"),
        ("file_signature", 'stat -c "%s|%Y|%y" -- "$1"',
         b"10|1729000000|2024-10-15 12:00:00.250000000 +0000\n",
         (10, 1729000000250000000)),
    ],
)
@pytest.mark.parametrize("returncode", (0, 3, 17))
def test_ssh_file_probes_use_same_absent_and_error_contract(
    method, command, stdout, expected, returncode,
):
    calls = []

    class ProbedTransport(SshTransport):
        async def run(self, argv, *, cwd=None, input=None, check=True):
            calls.append((argv, check))
            return subprocess.CompletedProcess(
                argv, returncode, stdout=stdout, stderr=b"probe failure",
            )

    path = "folder name/result.out"
    query = getattr(ProbedTransport("example-host"), method)
    if returncode == 17:
        with pytest.raises(subprocess.CalledProcessError) as error:
            asyncio.run(query(path))
        assert error.value.returncode == 17
        assert error.value.output == stdout
        assert error.value.stderr == b"probe failure"
    else:
        assert asyncio.run(query(path)) == (None if returncode == 3 else expected)
    assert calls == [
        (("sh", "-c", f'[ -f "$1" ] || exit 3; {command}', "sh", path), False),
    ]


@pytest.mark.parametrize("zstd", (False, True))
def test_ssh_batch_tar_uses_explicit_compression_without_probing(tmp_path, zstd):
    source = tmp_path / "source"
    source.mkdir()
    (source / "one.dat").write_bytes(b"content")
    destination = tmp_path / "remote"
    returned = tmp_path / "returned"
    transport = _LoopbackSshTransport(zstd=zstd)
    commands = []

    async def record(producer, consumer):
        commands.append((producer, consumer))

    transport._transfer_tar = record

    async def check():
        await transport.put_files("one.dat", source, destination)
        await transport.get_files("one.dat", source, returned)

    asyncio.run(check())
    assert len(commands) == 2
    for producer, consumer in commands:
        assert ("--zstd" in producer) is zstd
        assert ("--zstd" in consumer) is zstd


def test_ssh_endpoint_compression_roundtrips_from_transport():
    from ..endpoint import endpoint_for_executor, create_executor_from_endpoint
    from ..remote_endpoint import RemoteEndpointExecutor

    for compressed in (False, True):
        endpoint = endpoint_for_executor(
            RemoteEndpointExecutor(SshTransport("cluster", zstd=compressed))
        )
        assert endpoint.ssh_zstd is compressed
        restored = create_executor_from_endpoint(endpoint)
        assert restored.transport.zstd is compressed
