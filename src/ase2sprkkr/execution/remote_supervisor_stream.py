"""Transient single-socket streaming owned by the direct process supervisor.

Files are spooled only when explicitly requested. A missing stream file means
LIVE-only observation, not an implicit stdout/stderr artifact. The target owns
no run history; terminal facts are sent on this socket and via a short-lived
acknowledgeable receipt.
"""
from __future__ import annotations

import asyncio
import json
import os
import shutil
import struct
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from ._durability import fsync_directory

MAGIC = b"A2SFOLL1"
_HEADER = struct.Struct("!cQQ")
_CODES = {"out": b"O", "err": b"E"}
_BACK = {code: kind for kind, code in _CODES.items()}
_CHUNK = 65536
_HIGH = 512 * 1024
_LOW = 128 * 1024
_MAX_REQUEST = 4096


def frame(kind, offset, data):
    return _HEADER.pack(_CODES[kind], offset, len(data)) + data


async def receive(reader):
    """Yield offset-bearing stream frames followed by terminal facts."""
    if await reader.readexactly(len(MAGIC)) != MAGIC:
        raise ValueError("Unknown supervisor stream protocol")
    while True:
        code, offset, length = _HEADER.unpack(await reader.readexactly(_HEADER.size))
        if length > (_MAX_REQUEST if code == b"T" else _CHUNK):
            raise ValueError("Oversized supervisor frame")
        payload = await reader.readexactly(length)
        if code == b"T":
            yield "terminal", 0, json.loads(payload)
            return
        if code not in _BACK:
            raise ValueError("Invalid supervisor channel")
        yield _BACK[code], offset, payload


@dataclass(eq=False)
class _Follower:
    writer: asyncio.StreamWriter
    offsets: dict
    replay: bool = True
    wake: asyncio.Event = field(default_factory=asyncio.Event)
    done: asyncio.Event = field(default_factory=asyncio.Event)


class SupervisorStreams:
    """One event loop writes the spools and fans out output to followers."""

    def __init__(self, stdout_path=None, stderr_path=None, *, detached=True,
                 socket_dir=None):
        self.paths = {"out": None if stdout_path is None else os.fspath(stdout_path),
                      "err": None if stderr_path is None else os.fspath(stderr_path)}
        self.merged = (stdout_path is not None and stderr_path is not None
                       and os.path.abspath(stdout_path) == os.path.abspath(stderr_path))
        if self.merged:
            del self.paths["err"]
        self.detached = detached
        self.connected = asyncio.Event()
        self.disconnected = asyncio.Event()
        self.cancelled = asyncio.Event()
        self._session_follower = None
        self._spools = {}
        self._sizes = {kind: 0 for kind in self.paths}
        self._followers = set()
        self._server = None
        self._terminal = None
        self._socket_dir = socket_dir
        self.socket_path = None

    async def start(self):
        for kind, path in self.paths.items():
            if path is not None:
                Path(path).parent.mkdir(parents=True, exist_ok=True)
                # Truncate the declared stream at launch. Keep the actual file
                # available during execution and after failure or cancellation.
                self._spools[kind] = open(path, "wb", buffering=0)
        if self._socket_dir is None:
            self._socket_dir = tempfile.mkdtemp(prefix="ase2sprkkr-follow-")
        self.socket_path = os.path.join(self._socket_dir, "socket")
        self._server = await asyncio.start_unix_server(self._connect, self.socket_path)
        os.chmod(self.socket_path, 0o600)
        return self.socket_path

    def publish(self, kind, data):
        if self._terminal is not None:
            raise RuntimeError("Output published after termination")
        if self.merged and kind == "err":
            kind = "out"
        if not data:
            return
        data = bytes(data)
        offset = self._sizes[kind]
        spool = self._spools.get(kind)
        if spool is not None:
            spool.write(data)
            spool.flush()
        self._sizes[kind] += len(data)
        wire = frame(kind, offset, data)
        for follower in tuple(self._followers):
            # A spooled channel can catch up from disk. A live-only channel
            # cannot: forward it even while another channel is replaying.
            if follower.replay and self.paths[kind] is not None:
                follower.wake.set()
                continue
            if follower.writer.transport.get_write_buffer_size() >= _HIGH:
                if self.paths[kind] is None:
                    follower.writer.transport.abort()
                else:
                    follower.replay = True
                    follower.wake.set()
                continue
            follower.writer.write(wire)
            follower.offsets[kind] += len(data)
            if follower.writer.transport.get_write_buffer_size() >= _HIGH:
                if self.paths[kind] is None:
                    # There is nowhere to replay bytes lost to backpressure.
                    follower.writer.transport.abort()
                else:
                    follower.replay = True
                    follower.wake.set()

    def sync(self):
        """Durability barrier before publishing successful terminal outcome."""
        for kind, spool in self._spools.items():
            spool.flush()
            os.fsync(spool.fileno())
            fsync_directory(Path(self.paths[kind]).parent)

    def terminal(self, outcome):
        self._terminal = dict(outcome)
        for follower in tuple(self._followers):
            # A client that stopped reading must not prevent the detached
            # supervisor from exiting. Already published bytes remain in the
            # named spool and can be recovered on a later attach.
            if follower.writer.transport.get_write_buffer_size() >= _HIGH:
                follower.writer.transport.abort()
            follower.wake.set()

    async def _connect(self, reader, writer):
        follower = None
        disconnect_watch = None
        try:
            raw = await reader.readline()
            if not raw or len(raw) > _MAX_REQUEST:
                raise ValueError("Malformed supervisor request")
            request = json.loads(raw)
            if request.get("command") == "cancel":
                self.cancelled.set()
                writer.write(b'{"accepted":true}\n')
                await writer.drain()
                return
            if request.get("command") != "follow":
                raise ValueError("Unknown supervisor command")
            offsets = request.get("offsets")
            if not isinstance(offsets, dict) or set(offsets) != set(self.paths):
                raise ValueError("Invalid stream offsets")
            for kind, offset in offsets.items():
                if offset is None:
                    # from_bytes=None means observe only newly arriving bytes.
                    offsets[kind] = self._sizes[kind]
                    continue
                if type(offset) is not int or offset < 0:
                    raise ValueError("Invalid stream offset")
                if self.paths[kind] is not None and offset > self._sizes[kind]:
                    raise ValueError("Offset is beyond the spool")
                if self.paths[kind] is None:
                    if offset != 0:
                        raise ValueError("Live-only output cannot be replayed")
                    offsets[kind] = self._sizes[kind]
            writer.transport.set_write_buffer_limits(high=_HIGH, low=_LOW)
            follower = _Follower(writer, offsets.copy())
            if not self.detached:
                if self._session_follower is not None:
                    raise ValueError("Attached execution already has a client")
                self._session_follower = follower
                self.connected.set()
                disconnect_watch = asyncio.create_task(reader.read(1))
                disconnect_watch.add_done_callback(
                    lambda _: (self.disconnected.set(), follower.wake.set())
                )
            self._followers.add(follower)
            writer.write(MAGIC)
            while True:
                if disconnect_watch is not None and disconnect_watch.done():
                    return
                if follower.replay:
                    await writer.drain()
                    for kind, path in self.paths.items():
                        if path is None:
                            if follower.offsets[kind] != self._sizes[kind]:
                                raise ConnectionError("Live-only follower missed output")
                            continue
                        offset = follower.offsets[kind]
                        if offset < self._sizes[kind]:
                            with open(path, "rb", buffering=0) as spool:
                                spool.seek(offset)
                                data = spool.read(min(_CHUNK, self._sizes[kind] - offset))
                            if not data:
                                raise RuntimeError("Supervisor spool truncated")
                            writer.write(frame(kind, offset, data))
                            follower.offsets[kind] += len(data)
                            if (self._terminal is not None
                                    and writer.transport.get_write_buffer_size() >= _HIGH):
                                writer.transport.abort()
                                return
                            await writer.drain()
                            break
                    else:
                        # No suspension between checking EOF and enabling LIVE.
                        follower.replay = False
                if self._terminal is not None and all(
                    follower.offsets[kind] == self._sizes[kind] for kind in self.paths
                ):
                    payload = json.dumps(self._terminal, separators=(",", ":")).encode()
                    writer.write(_HEADER.pack(b"T", 0, len(payload)) + payload)
                    await writer.drain()
                    return
                if not follower.replay:
                    await follower.wake.wait()
                    follower.wake.clear()
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError, asyncio.IncompleteReadError):
            pass
        finally:
            if follower is not None:
                self._followers.discard(follower)
                if follower is self._session_follower and self._terminal is None:
                    self.disconnected.set()
                follower.done.set()
            if disconnect_watch is not None:
                disconnect_watch.cancel()
            writer.close()
            try:
                await writer.wait_closed()
            except (BrokenPipeError, ConnectionResetError):
                pass

    async def finish(self):
        """Let active clients drain the replay before removing the socket."""
        if self._server is not None:
            self._server.close()
            await self._server.wait_closed()
        if self._followers:
            await asyncio.gather(*(f.done.wait() for f in tuple(self._followers)))
        for spool in self._spools.values():
            spool.close()
        if self._socket_dir is not None:
            shutil.rmtree(self._socket_dir, ignore_errors=True)


async def follow_socket(path, offsets, on_chunk):
    reader, writer = await asyncio.open_unix_connection(path)
    try:
        writer.write(json.dumps({"command": "follow", "offsets": offsets}).encode() + b"\n")
        await writer.drain()
        async for kind, offset, data in receive(reader):
            if kind == "terminal":
                return data
            on_chunk(kind, offset, data)
        raise ConnectionError("Supervisor ended without terminal outcome")
    finally:
        writer.close()
        await writer.wait_closed()
