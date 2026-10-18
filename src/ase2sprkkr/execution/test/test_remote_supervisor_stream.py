"""Regression tests for the transient supervisor FOLLOW service."""
import asyncio

import pytest

from ..remote_supervisor_stream import SupervisorStreams, follow_socket, frame


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
