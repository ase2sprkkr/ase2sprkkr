"""End-to-end supervisor lifecycle tests, independent of the SSH transport."""

import asyncio
import os
import sys

from .. import remote_process_worker
from ..remote_supervisor_stream import follow_socket


def _spec(tmp_path, socket_fd, *, detached, argv=None):
    return {
        "argv": argv or [sys.executable, "-c", "import time; time.sleep(20)"],
        "directory": str(tmp_path),
        "output_file": str(tmp_path / "stdout"),
        "error_file": str(tmp_path / "stderr"),
        "receipt_id": "a" * 32,
        "socket_fd": socket_fd,
        "detached": detached,
    }


def test_attached_run_stops_when_follow_connection_drops(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args: captured.append(args))
        rd, wr = os.pipe()
        task = asyncio.create_task(remote_process_worker._run(_spec(tmp_path, wr, detached=False)))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
        await writer.drain()
        assert await reader.readexactly(8) == b"A2SFOLL1"
        writer.close()
        await writer.wait_closed()
        assert await asyncio.wait_for(task, timeout=5) != 0
        assert captured[0][1].value == "cancelled"
    asyncio.run(check())


def test_detached_run_survives_follow_disconnect(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args: captured.append(args))
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            sys.executable, "-c",
            "import time; print('first', flush=True); time.sleep(.2); print('last', flush=True)",
        ])
        task = asyncio.create_task(remote_process_worker._run(spec))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
        await writer.drain()
        assert await reader.readexactly(8) == b"A2SFOLL1"
        writer.close()
        await writer.wait_closed()
        assert await asyncio.wait_for(task, timeout=5) == 0
        assert (tmp_path / "stdout").read_bytes() == b"first\nlast\n"
        assert captured[0][1].value == "finished"
    asyncio.run(check())


def test_cancel_is_delivered_to_supervisor_not_just_supervisor_pid(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args, **kwargs: captured.append(args))
        rd, wr = os.pipe()
        task = asyncio.create_task(remote_process_worker._run(
            _spec(tmp_path, wr, detached=True)
        ))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"cancel"}\n')
        await writer.drain()
        assert b'"accepted":true' in await reader.readline()
        writer.close()
        await writer.wait_closed()
        await asyncio.wait_for(task, 5)
        assert captured[0][1].value == "cancelled"
    asyncio.run(check())


def test_failed_detached_launch_clears_previous_output(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(
            remote_process_worker, "publish_outcome",
            lambda *args, **kwargs: captured.append(args),
        )
        (tmp_path / "stdout").write_bytes(b"previous")
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            str(tmp_path / "missing-command"),
        ])
        result = await remote_process_worker._run(spec)
        with os.fdopen(rd, "rb") as pipe:
            assert pipe.readline()
        assert result != 0
        assert (tmp_path / "stdout").read_bytes() == b""
        assert captured[0][1].value == "failed"
        assert not list(tmp_path.glob(".*ase2sprkkr-*"))
    asyncio.run(check())
