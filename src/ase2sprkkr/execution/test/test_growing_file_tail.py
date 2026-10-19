"""Tests for the shared growing-file observation contract."""
import asyncio
import os

from ..growing_file_tail import GrowingFileTail, follow_growing_files


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
