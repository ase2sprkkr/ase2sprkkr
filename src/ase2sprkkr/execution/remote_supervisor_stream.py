"""Transient, backpressure-aware stream service owned by a direct-process supervisor.

A single Unix socket belongs to the supervisor, not to the client or a
separate lifecycle agent. The supervisor is the only writer of the spools.
FOLLOW offsets refer to the first *not yet queued* byte of each stream.
The connection handler never polls spool files for changes.
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

MAGIC = b"A2SFOLL1"
_HEADER = struct.Struct("!cQQ")  # channel, absolute byte offset, payload size
_CODES = {"out": b"O", "err": b"E"}
_BACK = {v: k for k, v in _CODES.items()}
_CHUNK = 65536
_HIGH = 512 * 1024
_LOW = 128 * 1024
_MAX_REQUEST = 4096


def frame(kind, offset, data):
    data = bytes(data)
    return _HEADER.pack(_CODES[kind], offset, len(data)) + data


async def receive(reader):
    """Yield verified, offset-bearing stream frames, followed by terminal JSON."""
    if await reader.readexactly(len(MAGIC)) != MAGIC:
        raise ValueError("Unknown supervisor stream protocol")
    while True:
        header = await reader.readexactly(_HEADER.size)
        channel, offset, length = _HEADER.unpack(header)
        if length > _CHUNK and channel != b"T":
            raise ValueError("Oversized supervisor frame")
        if channel == b"T":
            if length > _MAX_REQUEST:
                raise ValueError("Oversized terminal fact")
            yield "terminal", 0, json.loads(await reader.readexactly(length))
            return
        if channel not in _BACK:
            raise ValueError("Unknown supervisor channel")
        yield _BACK[channel], offset, await reader.readexactly(length)


@dataclass(eq=False)
class _Follower:
    writer: asyncio.StreamWriter
    offsets: dict
    replay: bool = True
    wake: asyncio.Event = field(default_factory=asyncio.Event)


class SupervisorStreams:
    """Single-event-loop spool fanout. Only the supervisor may call publish()."""

    def __init__(self, stdout_path, stderr_path, *, detached=True):
        self.paths = {"out": os.fspath(stdout_path)}
        self.merged = os.path.abspath(stdout_path) == os.path.abspath(stderr_path)
        if not self.merged:
            self.paths["err"] = os.fspath(stderr_path)
        self.detached = detached
        self.connected = asyncio.Event()
        self.disconnected = asyncio.Event()
        self._session_follower = None
        self._spools = {}
        self._sizes = {kind: 0 for kind in self.paths}
        self._followers = set()
        self._server = None
        self._terminal = None
        self._socket_dir = None
        self.socket_path = None

    async def start(self):
        for kind, path in self.paths.items():
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            self._spools[kind] = open(path, "wb", buffering=0)
        # A private, disposable path: no target-side lifecycle database.
        self._socket_dir = tempfile.mkdtemp(prefix="ase2sprkkr-follow-")
        self.socket_path = os.path.join(self._socket_dir, "socket")
        self._server = await asyncio.start_unix_server(self._connect, self.socket_path)
        os.chmod(self.socket_path, 0o600)
        return self.socket_path

    def publish(self, kind, data):
        """Spool/flush/fanout is one non-yielding event-loop critical section."""
        if self._terminal is not None:
            raise RuntimeError("Output published after terminal state")
        if self.merged and kind == "err":
            kind = "out"
        data = bytes(data)
        if not data:
            return
        stream = self._spools[kind]
        offset = self._sizes[kind]
        stream.write(data)
        self._sizes[kind] += len(data)
        # Flush before making the new high-water mark visible to replay.
        stream.flush()
        wire = frame(kind, offset, data)
        for follower in tuple(self._followers):
            if follower.replay:
                follower.wake.set()
                continue
            # A LIVE follower has exactly one writer. If slow, stop enqueuing
            # live bytes; REPLAY resumes from its first unqueued offset.
            if follower.writer.transport.get_write_buffer_size() >= _HIGH:
                follower.replay = True
                follower.wake.set()
                continue
            follower.writer.write(wire)
            follower.offsets[kind] += len(data)
            if follower.writer.transport.get_write_buffer_size() >= _HIGH:
                follower.replay = True
                follower.wake.set()

    def terminal(self, outcome):
        """Called only after waitpid, both pipe pumps, and outcome publication."""
        if self._terminal is not None:
            raise RuntimeError("Terminal outcome already published")
        self._terminal = dict(outcome)
        for follower in tuple(self._followers):
            follower.wake.set()

    async def _connect(self, reader, writer):
        follower = None
        disconnect_watch = None
        try:
            request = await reader.readline()
            if not request or len(request) > _MAX_REQUEST:
                raise ValueError("Invalid FOLLOW request")
            data = json.loads(request)
            if data.get("command") != "follow":
                raise ValueError("Unknown supervisor request")
            offsets = data.get("offsets")
            if not isinstance(offsets, dict) or set(offsets) != set(self.paths):
                raise ValueError("Invalid stream offsets")
            if any(
                type(offsets[kind]) is not int or not 0 <= offsets[kind] <= self._sizes[kind]
                for kind in self.paths
            ):
                raise ValueError("Stream offset outside the current spool")
            writer.transport.set_write_buffer_limits(high=_HIGH, low=_LOW)
            # Register before yielding: no output between snapshot and replay.
            follower = _Follower(writer, offsets.copy())
            if not self.detached:
                if self._session_follower is not None:
                    raise ValueError("Non-detached run already has its session")
                self._session_follower = follower
                self.connected.set()
                # The peer sends nothing after its FOLLOW request. Detect EOF
                # even while the process produces no output.
                disconnect_watch = asyncio.create_task(reader.read(1))
                disconnect_watch.add_done_callback(lambda _: self.disconnected.set())
            self._followers.add(follower)
            writer.write(MAGIC)
            while True:
                if follower.replay:
                    # Drain prior queued LIVE data before replaying the rest.
                    await writer.drain()
                    for kind in self.paths:
                        offset = follower.offsets[kind]
                        if offset < self._sizes[kind]:
                            with open(self.paths[kind], "rb", buffering=0) as spool:
                                spool.seek(offset)
                                data = spool.read(min(_CHUNK, self._sizes[kind] - offset))
                            if not data:
                                raise RuntimeError("Spool was truncated during FOLLOW")
                            writer.write(frame(kind, offset, data))
                            follower.offsets[kind] += len(data)
                            # Yield to the pumps and apply socket backpressure.
                            await writer.drain()
                            break
                    else:
                        # No await between checking spool EOF and registering
                        # LIVE. An arriving publish therefore cannot be lost.
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
        except (BrokenPipeError, ConnectionResetError, asyncio.IncompleteReadError):
            pass
        finally:
            if follower is not None:
                self._followers.discard(follower)
                if follower is self._session_follower and self._terminal is None:
                    self.disconnected.set()
            if disconnect_watch is not None:
                disconnect_watch.cancel()
            writer.close()
            try:
                await writer.wait_closed()
            except (BrokenPipeError, ConnectionResetError):
                pass

    async def finish(self, *, grace=2.0):
        """Allow terminal frames to drain, but never let a client hold the job."""
        if self._server is not None:
            self._server.close()
            await self._server.wait_closed()
        if self._followers:
            try:
                await asyncio.wait_for(
                    asyncio.gather(*(f.writer.wait_closed() for f in tuple(self._followers))),
                    timeout=grace,
                )
            except (asyncio.TimeoutError, ConnectionError):
                pass
        for follower in tuple(self._followers):
            follower.writer.close()
        for spool in self._spools.values():
            spool.close()
        if self._socket_dir is not None:
            shutil.rmtree(self._socket_dir, ignore_errors=True)


async def follow_socket(path, offsets, on_chunk):
    """Follow one supervisor connection; on_chunk(kind, offset, bytes) is sync."""
    reader, writer = await asyncio.open_unix_connection(path)
    try:
        writer.write(json.dumps({"command": "follow", "offsets": offsets}).encode() + b"\n")
        await writer.drain()
        outcome = None
        async for kind, offset, data in receive(reader):
            if kind == "terminal":
                outcome = data
                break
            on_chunk(kind, offset, data)
        return outcome
    finally:
        writer.close()
        await writer.wait_closed()
