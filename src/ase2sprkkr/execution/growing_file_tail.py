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

    def read_available(self, paths=None):
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
                    data = source.read(self.chunk_size)
            except FileNotFoundError:
                continue
            self.offsets[kind] = offset + len(data)
            if data:
                chunks.append((kind, data))
        return chunks

    async def wait_for_change(self):
        """Wake on a file event or periodically reconcile missed notifications."""
        try:
            await asyncio.wait_for(self._changed.wait(), self.poll_interval)
        except asyncio.TimeoutError:
            pass
        self._changed.clear()
