"""Tests for the shared growing-file observation contract."""
import asyncio
import os
from threading import Event

import pytest

from ..growing_file_tail import (
    GrowingFileTail, follow_growing_files, follow_job_files,
)


def test_append_and_eof(tmp_path):
    async def check():
        path = tmp_path / "out"
        path.write_bytes(b"first")
        async with GrowingFileTail({"out": path}) as tail:
            assert tail.read_available() == [("out", b"first")]
            assert tail.read_available() == []
            path.open("ab").write(b"second")
            await tail.wait_for_change()
            assert tail.read_available() == [("out", b"second")]
            assert tail.read_available() == []
    asyncio.run(check())


def test_missing_then_created(tmp_path):
    async def check():
        path = tmp_path / "new"
        async with GrowingFileTail({"out": path}) as tail:
            assert tail.read_available() == []
            path.write_bytes(b"ready")
            await tail.wait_for_change()
            assert tail.read_available() == [("out", b"ready")]
    asyncio.run(check())


def test_rename_keeps_offset(tmp_path):
    async def check():
        temporary = tmp_path / "launching"
        final = tmp_path / "out"
        temporary.write_bytes(b"abc")
        async with GrowingFileTail({"out": temporary}, watch_paths=(final,)) as tail:
            assert tail.read_available() == [("out", b"abc")]
            os.replace(temporary, final)
            with final.open("ab") as destination:
                destination.write(b"def")
            assert tail.read_available({"out": final}) == [("out", b"def")]
    asyncio.run(check())


def test_live_only_starts_at_eof(tmp_path):
    async def check():
        path = tmp_path / "out"
        path.write_bytes(b"previous")
        async with GrowingFileTail({"out": path}, offsets={"out": None}) as tail:
            assert tail.read_available() == []
            with path.open("ab") as destination:
                destination.write(b"new")
            assert tail.read_available() == [("out", b"new")]
    asyncio.run(check())


def test_bounded_chunks(tmp_path):
    async def check():
        path = tmp_path / "out"
        path.write_bytes(b"123456")
        async with GrowingFileTail({"out": path}, chunk_size=2) as tail:
            assert tail.read_available() == [("out", b"12")]
            assert tail.read_available() == [("out", b"34")]
            assert tail.read_available() == [("out", b"56")]
            assert tail.read_available() == []
    asyncio.run(check())


def test_replaced_file_restarts_at_zero(tmp_path):
    async def check():
        path = tmp_path / "out"
        path.write_bytes(b"old")
        async with GrowingFileTail({"out": path}) as tail:
            assert tail.read_available() == [("out", b"old")]
            replacement = tmp_path / "replacement"
            replacement.write_bytes(b"new content")
            os.replace(replacement, path)
            assert tail.read_available() == [("out", b"new content")]
    asyncio.run(check())


def test_invalid_fallback_interval():
    try:
        GrowingFileTail({"out": "/tmp/example"}, poll_interval=0)
    except ValueError:
        pass
    else:
        raise AssertionError("Nonpositive poll intervals must be rejected")


def test_follow_growing_files_drains_both_channels_before_terminal(tmp_path):
    async def check():
        out, err = tmp_path / "out", tmp_path / "err"
        out.write_bytes(b"123456")
        err.write_bytes(b"xyz")
        observed = []
        checks = []

        async def terminal():
            checks.append(True)
            return "done"

        result = await follow_growing_files(
            {"out": out, "err": err},
            lambda kind, data: observed.append((kind, data)),
            terminal,
        )
        assert result == "done"
        assert observed == [("out", b"123456"), ("err", b"xyz")]
        assert len(checks) == 1

    asyncio.run(check())


def test_follow_growing_files_skips_undefined_channels(tmp_path):
    async def check():
        seen = []
        async def terminal():
            return "done"
        result = await follow_growing_files(
            {"out": None, "err": None},
            lambda kind, data: seen.append((kind, data)),
            terminal,
        )
        assert result == "done"
        assert seen == []
        assert list(tmp_path.iterdir()) == []

    asyncio.run(check())


def test_follow_job_files_drains_output_created_at_completion(tmp_path):
    async def check():
        path = tmp_path / "later.out"
        entered, release = Event(), Event()
        seen = []
        calls = []

        def wait():
            calls.append("wait")
            entered.set()
            assert release.wait(2)
            path.write_bytes(b"final bytes")
            return "done"

        following = asyncio.create_task(follow_job_files(
            {"out": path}, lambda kind, data: seen.append((kind, data)),
            wait, poll_interval=0.01,
        ))
        try:
            assert await asyncio.to_thread(entered.wait, 2)
        finally:
            release.set()
        assert await asyncio.wait_for(following, 3) == "done"
        assert calls == ["wait"]
        assert seen == [("out", b"final bytes")]

    asyncio.run(check())


def test_follow_job_files_maps_terminal_status(tmp_path):
    async def check():
        observed = []

        async def terminal(native_status):
            observed.append(native_status)
            return {"state": "finished", "return_code": native_status}

        return_value = await follow_job_files(
            {"out": tmp_path / "missing"}, lambda *_: None,
            lambda: 42, terminal=terminal, poll_interval=0.01,
        )
        assert return_value == {"state": "finished", "return_code": 42}
        assert observed == [42]

    asyncio.run(check())


def test_follow_job_files_propagates_wait_errors(tmp_path):
    def wait():
        raise RuntimeError("scheduler wait failed")

    async def check():
        with pytest.raises(RuntimeError, match="scheduler wait failed"):
            await asyncio.wait_for(
                follow_job_files(
                    {"out": tmp_path / "missing"}, lambda *_: None,
                    wait, poll_interval=0.01,
                ),
                3,
            )

    asyncio.run(check())


def test_follow_job_files_cancelled_observation_does_not_cancel_job(tmp_path):
    async def check():
        entered, release, exited = Event(), Event(), Event()

        def wait():
            entered.set()
            assert release.wait(2)
            exited.set()
            return "finished"

        follower = asyncio.create_task(follow_job_files(
            {"out": tmp_path / "missing"}, lambda *_: None,
            wait, poll_interval=0.01,
        ))
        try:
            assert await asyncio.to_thread(entered.wait, 2)
            follower.cancel()
            with pytest.raises(asyncio.CancelledError):
                await follower
            assert not exited.is_set()
        finally:
            release.set()
        assert await asyncio.to_thread(exited.wait, 2)

    asyncio.run(check())


def test_terminal_output_drain_is_bounded_when_producer_keeps_writing(tmp_path):
    async def check():
        path = tmp_path / "stream"
        path.write_bytes(b"A" * (65536 * 2))
        received = []
        checks = []
        def on_chunk(kind, data):
            received.append(data)
            # An ongoing producer appends more data for every chunk we read.
            with path.open("ab") as sink:
                sink.write(b"more data")

        async def terminal():
            checks.append(True)
            return "unknown"

        result = await asyncio.wait_for(follow_growing_files(
            {"out": path}, on_chunk, terminal, poll_interval=0.01,
        ), 3)
        assert result == "unknown"
        assert 2 <= len(received) <= 3
        assert len(checks) == 1

    asyncio.run(check())



def test_read_error_does_not_wait_for_scheduler_job_to_end(tmp_path, monkeypatch):
    from threading import Event

    entered, release = Event(), Event()

    def wait():
        entered.set()
        release.wait(10)
        return "finished"

    def unreadable(self, *args, **kwargs):
        raise PermissionError("output inaccessible")

    monkeypatch.setattr(GrowingFileTail, "read_available", unreadable)

    async def check():
        with pytest.raises(PermissionError, match="inaccessible"):
            await asyncio.wait_for(
                follow_job_files(
                    {"out": tmp_path / "job.out"}, lambda *_: None,
                    wait, poll_interval=0.01,
                ), 2,
            )
        assert await asyncio.to_thread(entered.wait, 1)
        assert not release.is_set()

    # Release only after asyncio.run has shut down its executor; this would
    # block with the former to_thread implementation.
    try:
        asyncio.run(check())
    finally:
        release.set()


def test_terminal_creates_output_after_empty_read(tmp_path):
    async def check():
        path = tmp_path / "late.out"
        seen = []

        async def terminal():
            path.write_bytes(b"last bytes")
            return "finished"

        result = await asyncio.wait_for(
            follow_growing_files(
                {"out": path}, lambda kind, data: seen.append((kind, data)),
                terminal, poll_interval=0.01,
            ), 2,
        )
        assert result == "finished"
        assert seen == [("out", b"last bytes")]

    asyncio.run(check())



@pytest.mark.parametrize("change", ("remove", "truncate"))
def test_terminal_drain_ends_if_snapshotted_bytes_disappear(
    tmp_path, monkeypatch, change,
):
    """Simulate a file disappearing *after* the terminal size was frozen."""
    async def check():
        path = tmp_path / "out"
        path.write_bytes(b"X" * (65536 * 3))
        original_snapshot = GrowingFileTail.snapshot_sizes
        snapshots = []

        def snapshot(tail):
            sizes = original_snapshot(tail)
            snapshots.append(sizes)
            if len(snapshots) == 1:
                # Only the initial terminal snapshot triggers the race.
                # A later size probe must observe the missing/truncated file.
                if change == "remove":
                    path.unlink()
                else:
                    path.write_bytes(b"")
            return sizes

        monkeypatch.setattr(GrowingFileTail, "snapshot_sizes", snapshot)

        async def terminal():
            return "finished"

        seen = []
        assert await asyncio.wait_for(follow_growing_files(
            {"out": path}, lambda kind, data: seen.append(data),
            terminal, poll_interval=0.01,
        ), 1) == "finished"
        assert snapshots[0] == {"out": 65536 * 3}
        assert snapshots[-1] == {"out": 0}
        assert seen == [b"X" * 65536]

    asyncio.run(check())


def test_continuous_output_yields_to_scheduler_completion(tmp_path):
    async def check():
        path = tmp_path / "busy"
        path.write_bytes(b"X" * (65536 * 10))
        completed = asyncio.Event()
        received = []

        def on_chunk(kind, data):
            received.append(data)
            with path.open("ab") as stream:
                stream.write(data)

        async def mark_terminal():
            await asyncio.sleep(0)
            completed.set()

        asyncio.create_task(mark_terminal())

        async def terminal():
            return "finished" if completed.is_set() else None

        result = await asyncio.wait_for(follow_growing_files(
            {"out": path}, on_chunk, terminal, poll_interval=0.01,
        ), 1)
        assert result == "finished"
        # Without an actual event-loop yield an ever-growing file can keep
        # terminal notification queued forever. Allow scheduling variations.
        assert len(received) < 100

    asyncio.run(check())


def test_callback_wakes_silent_file_follower_without_poll_timeout(tmp_path):
    async def check():
        loop = asyncio.get_running_loop()
        completed = loop.create_future()
        entered = asyncio.Event()
        seen = []

        async def terminal():
            entered.set()
            return "finished" if completed.done() else None

        task = asyncio.create_task(follow_growing_files(
            {"out": tmp_path / "does-not-exist"}, lambda *args: seen.append(args),
            terminal, poll_interval=30, wakeup=completed,
        ))
        await asyncio.wait_for(entered.wait(), 1)
        completed.set_result(True)
        assert await asyncio.wait_for(task, 2) == "finished"
        assert seen == []

    asyncio.run(check())
