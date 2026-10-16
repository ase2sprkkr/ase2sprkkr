"""Framing for calculation stdout/stderr over one transport data channel.

The protocol preserves stdout/stderr identity while the streams remain
separate. If a caller explicitly requests a merged sink (``error_file=True``
at the SPR-KKR API), the merge is intentionally lossy: all data is retained in
the output, but cross-stream ordering is not guaranteed to be reproducible and
the original channel identity cannot be recovered afterwards.
"""

from __future__ import annotations

import asyncio
import struct


STREAM_MAGIC = b"A2SSTRM1"
_HEADER = struct.Struct("!cQ")
_KIND_TO_CODE = {"out": b"O", "err": b"E"}
_CODE_TO_KIND = {value: key for key, value in _KIND_TO_CODE.items()}


class StreamProtocolError(RuntimeError):
    """Malformed or unsupported remote stream protocol data."""


def encode_stream_frame(kind, data):
    try:
        code = _KIND_TO_CODE[kind]
    except KeyError as exc:
        raise ValueError(f"Unknown stream kind: {kind!r}") from exc
    data = bytes(data)
    return _HEADER.pack(code, len(data)) + data


async def _read_exactly(source, size):
    try:
        return await source.readexactly(size)
    except asyncio.IncompleteReadError as exc:
        raise StreamProtocolError(
            "Remote stream ended in the middle of a frame"
        ) from exc


async def read_stream_frames(source):
    """Yield ``(kind, data)`` records from a framed remote stream."""
    magic = await _read_exactly(source, len(STREAM_MAGIC))
    if magic != STREAM_MAGIC:
        raise StreamProtocolError("Remote stream has an unknown protocol header")

    while True:
        first = await source.read(1)
        if not first:
            return
        rest = await _read_exactly(source, _HEADER.size - 1)
        code, size = _HEADER.unpack(first + rest)
        try:
            kind = _CODE_TO_KIND[code]
        except KeyError as exc:
            raise StreamProtocolError(
                f"Remote stream contains an unknown channel code: {code!r}"
            ) from exc
        yield kind, await _read_exactly(source, size)


class FramedRunningCommand:
    """Expose calculation streams from one framed transport command.

    The wrapped process' stdout is the framed data channel. Its stderr is a
    transport/control diagnostic channel and is deliberately never exposed as
    calculation stderr.
    """

    def __init__(self, process, *, diagnostic_limit=65536):
        self._process = process
        self.stdout = asyncio.StreamReader()
        self.stderr = asyncio.StreamReader()
        self.returncode = None
        self._diagnostic_limit = diagnostic_limit
        self._diagnostics = bytearray()
        self._protocol_error = None
        self._frames_task = asyncio.create_task(self._demultiplex())
        self._diagnostics_task = asyncio.create_task(self._drain_diagnostics())

    @property
    def diagnostics(self):
        return bytes(self._diagnostics)

    async def _demultiplex(self):
        try:
            async for kind, data in read_stream_frames(self._process.stdout):
                (self.stderr if kind == "err" else self.stdout).feed_data(data)
        except asyncio.CancelledError:
            raise
        except Exception as exc:
            self._protocol_error = exc
            try:
                self._process.terminate()
            except (ProcessLookupError, RuntimeError):
                pass
        finally:
            self.stdout.feed_eof()
            self.stderr.feed_eof()

    async def _drain_diagnostics(self):
        while True:
            data = await self._process.stderr.read(65536)
            if not data:
                return
            remaining = self._diagnostic_limit - len(self._diagnostics)
            if remaining > 0:
                self._diagnostics.extend(data[:remaining])

    async def wait(self):
        self.returncode = await self._process.wait()
        await asyncio.gather(self._frames_task, self._diagnostics_task)
        diagnostic = self.diagnostics.decode("utf8", errors="replace").strip()
        if self._protocol_error is not None:
            message = str(self._protocol_error)
            if diagnostic:
                message += f"; remote transport diagnostic: {diagnostic}"
            raise StreamProtocolError(message) from self._protocol_error
        if self.returncode != 0:
            message = f"Remote stream helper exited with status {self.returncode}"
            if diagnostic:
                message += f"; remote transport diagnostic: {diagnostic}"
            raise StreamProtocolError(message)
        return self.returncode

    def terminate(self):
        return self._process.terminate()
