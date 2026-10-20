"""Incremental observation of growing files, shared by local and remote runners.

Like GNU tail -f, use inotify where available and recheck periodically where
filesystem notifications are unsupported or not forwarded (notably some network
filesystems). Ordinary file EOF is not a waitable event, even under asyncio.
"""
from __future__ import annotations

import asyncio
import ctypes
import os
import sys
from pathlib import Path
from threading import Thread


# Directory watches also detect an output file created after observation starts
# or atomically renamed into place.
_WATCH_MASK = (0x00000002 | 0x00000004 | 0x00000008 | 0x00000080 |
               0x00000100 | 0x00000200 | 0x00000400 | 0x00000800)


class GrowingFileTail:
    """Read only new bytes, and await changes without busy-polling.

    Offsets are per logical channel, not pathname: a worker can atomically
    rename its launching spool to the final output without replaying it.
    None as an initial offset means start at current EOF (live-only observation).
    """

    def __init__(self, paths, *, offsets=None, watch_paths=(),
                 poll_interval=1.0, chunk_size=65536):
        if poll_interval <= 0:
            raise ValueError("poll_interval must be positive")
        self.paths = dict(paths)
        self.offsets = offsets if offsets is not None else {
            kind: 0 for kind in self.paths
        }
        self.watch_paths = tuple(watch_paths)
        self.poll_interval = poll_interval
        self.chunk_size = chunk_size
        self._identities = {}
        self._fd = None
        self._loop = None
        self._changed = asyncio.Event()

    async def __aenter__(self):
        self._loop = asyncio.get_running_loop()
        if sys.platform != "linux" or not hasattr(self._loop, "add_reader"):
            return self
        libc = ctypes.CDLL(None, use_errno=True)
        init = getattr(libc, "inotify_init1", None)
        add_watch = getattr(libc, "inotify_add_watch", None)
        if init is None or add_watch is None:
            return self
        init.argtypes = [ctypes.c_int]
        init.restype = ctypes.c_int
        add_watch.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
        add_watch.restype = ctypes.c_int
        fd = init(os.O_NONBLOCK | os.O_CLOEXEC)
        if fd < 0:
            return self
        try:
            watched = set()
            for path in (*self.paths.values(), *self.watch_paths):
                if path is None:
                    continue
                directory = os.path.abspath(os.fspath(Path(path).parent))
                if directory in watched:
                    continue
                watched.add(directory)
                add_watch(fd, os.fsencode(directory), _WATCH_MASK)
            # If no directories can be watched, rely on the fallback timer.
            self._loop.add_reader(fd, self._notify)
        except (OSError, NotImplementedError):
            os.close(fd)
        else:
            self._fd = fd
        return self

    async def __aexit__(self, *_):
        if self._fd is not None:
            self._loop.remove_reader(self._fd)
            os.close(self._fd)
            self._fd = None

    def _notify(self):
        try:
            while True:
                if not os.read(self._fd, 65536):
                    break
        except BlockingIOError:
            pass
        self._changed.set()

    def snapshot_sizes(self):
        """Freeze bytes already present when job completion is observed."""
        sizes = {}
        for kind, path in self.paths.items():
            try:
                sizes[kind] = 0 if path is None else os.path.getsize(path)
            except FileNotFoundError:
                sizes[kind] = 0
        return sizes

    def read_available(self, paths=None, *, max_sizes=None):
        """Return bounded (channel, bytes) chunks; never block at EOF."""
        paths = self.paths if paths is None else paths
        chunks = []
        seen = set()
        for kind, path in paths.items():
            if path is None:
                continue
            path = os.fspath(path)
            if path in seen:
                continue
            seen.add(path)
            try:
                with open(path, "rb", buffering=0) as source:
                    stat = os.fstat(source.fileno())
                    size = stat.st_size
                    identity = (stat.st_dev, stat.st_ino)
                    previous = self._identities.get(kind)
                    offset = self.offsets[kind]
                    if previous is not None and previous != identity:
                        offset = 0
                    self._identities[kind] = identity
                    if offset is None:
                        self.offsets[kind] = size
                        continue
                    if size < offset:  # truncated or replaced in place
                        offset = 0
                    source.seek(offset)
                    remaining = (
                        self.chunk_size if max_sizes is None else
                        min(self.chunk_size, max(0, max_sizes[kind] - offset))
                    )
                    data = source.read(remaining)
            except FileNotFoundError:
                continue
            self.offsets[kind] = offset + len(data)
            if data:
                chunks.append((kind, data))
        return chunks

    async def wait_for_change(self, *, wakeup=None):
        """Wake on file changes, optional terminal future, or filesystem timeout."""
        if wakeup is None:
            try:
                await asyncio.wait_for(self._changed.wait(), self.poll_interval)
            except asyncio.TimeoutError:
                pass
        else:
            file_change = asyncio.create_task(self._changed.wait())
            try:
                await asyncio.wait(
                    (file_change, wakeup), timeout=self.poll_interval,
                    return_when=asyncio.FIRST_COMPLETED,
                )
            finally:
                file_change.cancel()
                await asyncio.gather(file_change, return_exceptions=True)
        self._changed.clear()

async def follow_growing_files(paths, on_chunk, terminal, *, offsets=None,
                               poll_interval=1.0, wakeup=None):
    """Observe file bytes and drain them before returning a terminal fact.

    `terminal` is an async callback returning None while execution remains
    active, or a terminal value when the producer has finished writing.
    This helper does not own backend status, files or transport.
    """
    async with GrowingFileTail(
        paths, offsets=offsets, poll_interval=poll_interval,
    ) as tail:
        finished = None
        limits = None
        while True:
            chunks = tail.read_available(max_sizes=limits)
            for kind, data in chunks:
                on_chunk(kind, data)
            if finished is None:
                # A producer may write indefinitely after contact is lost.
                # Check the wait outcome after every bounded read.
                finished = await terminal()
                if finished is not None:
                    limits = tail.snapshot_sizes()
            if finished is not None:
                # Drain bytes captured at terminal observation, including a
                # file that appeared *between* the read and size snapshot.
                pending = [
                    (kind, size) for kind, size in limits.items()
                    if size and (
                        tail.offsets[kind] is None or tail.offsets[kind] < size
                    )
                ]
                if not pending:
                    return finished
                if not chunks:
                    # A snapshotted file may subsequently disappear/truncate.
                    # In that case its missing bytes cannot be drained. Do not
                    # spin on the frozen size, but do retry files still there.
                    available = tail.snapshot_sizes()
                    if all(available[kind] < size for kind, size in pending):
                        return finished
            if chunks:
                # read_available(), terminal() and on_chunk() may all finish
                # synchronously: explicitly let scheduler callbacks run.
                await asyncio.sleep(0)
            else:
                await tail.wait_for_change(wakeup=wakeup)

async def follow_job_files(paths, on_chunk, wait, *, terminal=None, offsets=None,
                           poll_interval=1.0):
    """Follow scheduler output while one blocking job.wait() runs in a thread.

    The optional async terminal callback maps the native job status to the
    caller's terminal facts. The common follower drains output before returning
    those facts. Cancelling observation does not cancel the scheduler job.
    """
    # A failed or cancelled observation must not make asyncio.run() wait
    # for a still-running scheduler calculation. Keep wait in a daemon thread.
    loop = asyncio.get_running_loop()
    waiting = loop.create_future()

    def deliver(value, error):
        if waiting.done():
            return
        if error is not None:
            waiting.set_exception(error)
        else:
            waiting.set_result(value)

    def invoke_wait():
        try:
            value, error = wait(), None
        except BaseException as exc:
            value, error = None, exc
        try:
            loop.call_soon_threadsafe(deliver, value, error)
        except RuntimeError:
            if not loop.is_closed():
                raise

    Thread(target=invoke_wait, daemon=True, name="psij-job-wait").start()

    async def outcome():
        if not waiting.done():
            return None
        status = await waiting
        return await terminal(status) if terminal is not None else status

    try:
        return await follow_growing_files(
            paths, on_chunk, outcome, offsets=offsets,
            poll_interval=poll_interval,
        )
    finally:
        if not waiting.done():
            waiting.cancel()
