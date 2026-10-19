"""Regression tests for the transient supervisor FOLLOW service."""
import asyncio

import pytest

from ..remote_supervisor_stream import SupervisorStreams, _Follower, _HIGH, follow_socket, frame


def test_supervisor_live_replay_and_terminal(tmp_path):
    async def run():
        output = tmp_path / "out"
        error = tmp_path / "err"
        streams = SupervisorStreams(output, error)
        await streams.start()
        try:
            streams.publish("out", b"before\n")
            streams.publish("err", b"warning\n")
            received = []

            async def follow():
                return await follow_socket(
                    streams.socket_path, {"out": 0, "err": 0},
                    lambda kind, offset, data: received.append((kind, offset, data)),
                )

            follower = asyncio.create_task(follow())
            # FOLLOW returns the backlog before transitioning to LIVE.
            await asyncio.sleep(0.05)
            streams.publish("out", b"after\n")
            streams.terminal({"state": "finished", "return_code": 0})
            outcome = await asyncio.wait_for(follower, 2)
            assert outcome == {"state": "finished", "return_code": 0}
            assert received == [
                ("out", 0, b"before\n"),
                ("err", 0, b"warning\n"),
                ("out", 7, b"after\n"),
            ]
            assert output.read_bytes() == b"before\nafter\n"
            assert error.read_bytes() == b"warning\n"
        finally:
            await streams.finish()
    asyncio.run(run())


def test_supervisor_follow_resumes_from_offset(tmp_path):
    async def run():
        streams = SupervisorStreams(tmp_path / "out", tmp_path / "err")
        await streams.start()
        try:
            streams.publish("out", b"abcdef")
            streams.publish("err", b"12345")
            got = []
            follower = asyncio.create_task(
                follow_socket(
                    streams.socket_path, {"out": 3, "err": 2},
                    lambda *args: got.append(args),
                )
            )
            await asyncio.sleep(0.05)
            streams.terminal({"state": "failed", "return_code": 7})
            result = await asyncio.wait_for(follower, 2)
            assert got == [("out", 3, b"def"), ("err", 2, b"345")]
            assert result["return_code"] == 7
        finally:
            await streams.finish()
    asyncio.run(run())


def test_supervisor_merged_output_uses_one_cursor(tmp_path):
    async def run():
        target = tmp_path / "merged"
        streams = SupervisorStreams(target, target)
        await streams.start()
        try:
            streams.publish("out", b"stdout")
            streams.publish("err", b"stderr")
            got = []
            follower = asyncio.create_task(
                follow_socket(
                    streams.socket_path, {"out": 0},
                    lambda *args: got.append(args),
                )
            )
            await asyncio.sleep(0.05)
            streams.terminal({"state": "finished", "return_code": 0})
            await asyncio.wait_for(follower, 2)
            assert got == [("out", 0, b"stdoutstderr")]
            assert target.read_bytes() == b"stdoutstderr"
        finally:
            await streams.finish()
    asyncio.run(run())


def test_supervisor_rejects_offsets_outside_spool(tmp_path):
    async def run():
        streams = SupervisorStreams(tmp_path / "out", tmp_path / "err")
        await streams.start()
        try:
            streams.publish("out", b"a")
            with pytest.raises((asyncio.IncompleteReadError, ConnectionError)):
                await follow_socket(
                    streams.socket_path,
                    {"out": 2, "err": 0},
                    lambda *_: None,
                )
        finally:
            await streams.finish()
    asyncio.run(run())


def test_attached_supervisor_detects_session_disconnect(tmp_path):
    async def run():
        streams = SupervisorStreams(tmp_path / "out", tmp_path / "err", detached=False)
        await streams.start()
        reader = writer = None
        try:
            reader, writer = await asyncio.open_unix_connection(streams.socket_path)
            writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
            await writer.drain()
            await asyncio.wait_for(streams.connected.wait(), 2)
            writer.close()
            await writer.wait_closed()
            await asyncio.wait_for(streams.disconnected.wait(), 2)
        finally:
            if writer is not None:
                writer.close()
            await streams.finish()
    asyncio.run(run())


def test_detached_supervisor_does_not_cancel_on_disconnect(tmp_path):
    async def run():
        streams = SupervisorStreams(tmp_path / "out", tmp_path / "err", detached=True)
        await streams.start()
        try:
            reader, writer = await asyncio.open_unix_connection(streams.socket_path)
            writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
            await writer.drain()
            writer.close()
            await writer.wait_closed()
            await asyncio.sleep(0.05)
            assert not streams.disconnected.is_set()
            streams.terminal({"state": "finished", "return_code": 0})
        finally:
            await streams.finish()
    asyncio.run(run())


def test_none_streams_live_only_without_spool_files(tmp_path):
    async def run():
        streams = SupervisorStreams(None, None)
        await streams.start()
        try:
            received = []
            task = asyncio.create_task(follow_socket(
                streams.socket_path, {"out": 0, "err": 0},
                lambda *args: received.append(args),
            ))
            await asyncio.sleep(0.05)
            streams.publish("out", b"observed")
            streams.publish("err", b"warning")
            streams.terminal({"state": "finished", "return_code": 0})
            assert (await asyncio.wait_for(task, 2))["return_code"] == 0
            assert [item[2] for item in received] == [b"observed", b"warning"]
            assert list(tmp_path.iterdir()) == []
            assert streams._spools == {}
        finally:
            await streams.finish()
    asyncio.run(run())


def test_live_only_stream_does_not_replay_earlier_bytes(tmp_path):
    async def run():
        streams = SupervisorStreams(None, None)
        await streams.start()
        try:
            streams.publish("out", b"old")
            seen = []
            follower = asyncio.create_task(follow_socket(
                streams.socket_path, {"out": 0, "err": 0},
                lambda *args: seen.append(args),
            ))
            await asyncio.sleep(0.05)
            streams.publish("out", b"new")
            streams.terminal({"state": "finished", "return_code": 0})
            await asyncio.wait_for(follower, 2)
            assert seen == [("out", 3, b"new")]
        finally:
            await streams.finish()
    asyncio.run(run())



def test_supervisor_writes_declared_output_during_execution(tmp_path):
    async def run():
        output, error = tmp_path / "out", tmp_path / "err"
        output.write_bytes(b"old")
        streams = SupervisorStreams(output, error)
        await streams.start()
        try:
            assert output.read_bytes() == b""
            streams.publish("out", b"partial")
            streams.publish("err", b"diagnostic")
            assert output.read_bytes() == b"partial"
            assert error.read_bytes() == b"diagnostic"
            streams.terminal({"state": "cancelled", "return_code": -15})
        finally:
            await streams.finish()
        assert output.read_bytes() == b"partial"
        assert error.read_bytes() == b"diagnostic"
    asyncio.run(run())


def test_no_stream_files_are_created(tmp_path):
    async def run():
        streams = SupervisorStreams(None, None)
        await streams.start()
        try:
            streams.publish("out", b"live only")
            streams.sync()
            assert streams._spools == {}
            assert list(tmp_path.iterdir()) == []
        finally:
            await streams.finish()
    asyncio.run(run())



def test_terminal_drops_backpressured_client_without_discarding_output(tmp_path):
    async def run():
        output = tmp_path / "out"
        streams = SupervisorStreams(output, None)
        await streams.start()
        writer = None
        try:
            reader, writer = await asyncio.open_unix_connection(streams.socket_path)
            writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
            await writer.drain()
            assert await reader.readexactly(8) == b"A2SFOLL1"
            for _ in range(1024):
                streams.publish("out", b"x" * 8192)
            streams.terminal({"state": "finished", "return_code": 0})
            await asyncio.wait_for(streams.finish(), 4)
            assert output.stat().st_size == 1024 * 8192
        finally:
            if writer is not None:
                writer.close()
    asyncio.run(run())


class _FakeFollowerTransport:
    def __init__(self, buffered=0):
        self.buffered = buffered
        self.aborted = False

    def get_write_buffer_size(self):
        return self.buffered

    def abort(self):
        self.aborted = True


class _FakeFollowerWriter:
    def __init__(self, buffered=0):
        self.transport = _FakeFollowerTransport(buffered)
        self.frames = []

    def write(self, data):
        self.frames.append(data)


def test_live_only_channel_is_delivered_during_other_channel_replay():
    async def run():
        # A follower catching up on a saved stdout must not lose live stderr.
        streams = SupervisorStreams("saved-stdout", None)
        streams.publish("out", b"previous stdout")
        writer = _FakeFollowerWriter()
        follower = _Follower(writer, {"out": 0, "err": 0}, replay=True)
        streams._followers.add(follower)

        streams.publish("err", b"live stderr")

        assert follower.replay
        assert writer.frames == [frame("err", 0, b"live stderr")]
        assert follower.offsets == {"out": 0, "err": len(b"live stderr")}
        assert not writer.transport.aborted

    asyncio.run(run())


def test_backpressured_live_only_channel_disconnects_instead_of_losing_bytes():
    async def run():
        streams = SupervisorStreams("saved-stdout", None)
        writer = _FakeFollowerWriter(buffered=_HIGH)
        follower = _Follower(writer, {"out": 0, "err": 0}, replay=True)
        streams._followers.add(follower)

        streams.publish("err", b"must not vanish")

        assert writer.transport.aborted
        assert follower.offsets["err"] == 0
        assert writer.frames == []

    asyncio.run(run())


def test_supervisor_finish_is_bounded_when_client_never_drains(tmp_path):
    async def run():
        streams = SupervisorStreams(tmp_path / "out", None)
        await streams.start()

        class Transport:
            def __init__(self, done):
                self.done = done
                self.aborted = False

            def get_write_buffer_size(self):
                return 0

            def abort(self):
                self.aborted = True
                self.done.set()

        class Writer:
            def __init__(self, transport):
                self.transport = transport

        done = asyncio.Event()
        transport = Transport(done)
        follower = _Follower(Writer(transport), {"out": 0, "err": 0})
        streams._followers.add(follower)
        streams.terminal({"state": "finished", "return_code": 0})
        try:
            await asyncio.wait_for(streams.finish(), 4)
            assert transport.aborted
            assert streams._spools["out"].closed
        finally:
            streams._followers.discard(follower)

    asyncio.run(run())
