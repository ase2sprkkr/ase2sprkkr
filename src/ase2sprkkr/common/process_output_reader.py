"""Parsers and stream helpers for calculation output."""

import asyncio
import functools
import os
import numpy as np
from .async_utils import run_coro_sync
from .decorators import maybeclassmethod
from .file_utils import filename_from_file
from .output_printing import normalize_print_output


class ProcessOutputParser:
    """Parse asynchronous stdout/stderr streams.

    Process creation, cancellation, file writing and raw-output delivery belong
    to the execution backend. Descendants only interpret the supplied streams.
    """

    def __init__(self, print_output=False):
        self.set_print_output(print_output)

    def result(self, output, error, return_code):
        """Postprocess parsed stdout/stderr after execution completes."""
        if return_code != 0:
            exc = ValueError(f"The process ended with return value {return_code}")
            exc.output = output
            exc.error = error
            raise exc
        return output, error, return_code

    async def read_error(self, stderr, *args):
        while True:
            line = await stderr.readline()
            if not line:
                return

    async def read_output(self, stdout, *args):
        raise NotImplementedError("Please, redefine ProcessOutputParser.read_output coroutine")

    def set_print_output(self, print_output):
        show_all, self.print_info, patterns = normalize_print_output(print_output)
        self.print_output = show_all or bool(patterns)

    @maybeclassmethod
    def read_from_file(self, cls, output, error=None, read_args=[], return_code=0):
        """Synchronous wrapper for reading output/error from files."""
        if not self:
            self = cls()
        return run_coro_sync(self.read_output_file(output, error, read_args, return_code))

    async def read_output_file(self, output, error=None, read_args=[], return_code=0):
        """Read output and error files through the same async parser API."""

        async def read_output_file():
            with AsyncioFileReader(output) as air:
                return await self.read_output(air, *read_args)

        async def read_error_file():
            if not error:
                return None
            with AsyncioFileReader(error) as air:
                return await self.read_error(air, *read_args)

        out_result, err_result = await asyncio.gather(read_output_file(), read_error_file())
        return self.result(out_result, err_result, return_code)


class AsyncioFileReader:
    """File reader that mimics asyncio StreamReader.

    It is synchronous internally, but exposes the async methods used by output
    parsers so the same parser can consume a subprocess stream or a saved file.
    """

    def __init__(self, filename, buffersize=8192):
        self._owns_file = isinstance(filename, (str, bytes, os.PathLike))
        self.file = None
        self.file = open(os.fspath(filename), "rb") if self._owns_file else filename
        self.buffersize = buffersize

    def __fspath__(self):
        filename = filename_from_file(self.file, None)
        if filename is None:
            raise TypeError("The underlying stream has no filesystem path")
        return filename

    def close(self):
        file, self.file = self.file, None
        if file is not None and self._owns_file:
            file.close()

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    async def __aiter__(self):
        while (line := self.file.readline()) != b"":
            yield line

    def __anext__(self):
        return self.readline()

    async def readline(self):
        return self.file.readline()

    @staticmethod
    @functools.lru_cache(maxsize=128)
    def separator_shifts(sep):
        shifts = np.ones(len(sep) + 1, dtype=int)
        shift = 1
        for pos in range(len(sep)):
            while shift <= pos and sep[pos] != sep[pos - shift]:
                shift += shifts[pos - shift]
            shifts[pos + 1] = shift
        return shifts

    async def readuntil(self, sep=b"\n"):
        lsep = len(sep)
        if lsep == 0:
            return b""

        shifts = self.separator_shifts(sep)

        out = b""
        buffer = bytearray(self.buffersize + lsep)
        bufEnd = 0
        bufPos = 0

        def result():
            if not out:
                return bytes(buffer[:bufPos])
            return out + bytes(buffer[:bufPos])

        startPos = 0
        matchLen = 0
        while True:
            if bufPos == bufEnd:
                if bufEnd >= self.buffersize:
                    out += buffer[:bufEnd]
                    bufEnd = bufPos = 0
                n = lsep - matchLen
                data = self.file.read(n)
                ldata = len(data)
                buffer[bufPos : bufPos + ldata] = data
                bufEnd += ldata
                if len(data) < n:
                    raise asyncio.IncompleteReadError(result(), sep)
            c = buffer[bufPos]
            bufPos += 1
            while matchLen >= 0 and sep[matchLen] != c:
                startPos += shifts[matchLen]
                matchLen -= shifts[matchLen]
            matchLen += 1
            if matchLen == lsep:
                return result()


async def readline(stdout):
    line = await stdout.readline()
    if not line:
        raise EOFError()
    return line.decode("utf8")


async def readline_until(stdout, cond, can_end=True):
    while True:
        line = await stdout.readline()
        if not line:
            if can_end:
                return None
            raise EOFError()
        if cond(line):
            return line.decode("utf8")


