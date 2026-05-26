"""Readers for processes outputs.

They employ asyncio to makes possible read stdio and stderr concurently.
"""

import asyncio
import functools
import subprocess
import os
import numpy as np
from threading import Thread
from .decorators import maybeclassmethod


class ProcessOutputParser:
    """
    Class, that run a process, optionally saves all the output of the process to a file,
    and pass the stdout and stderr of the process to its two async routines, read_error
    and read_output.

    The descendant can redefine the routines to parse the output (or its parts).
    """

    def __init__(self, print_output=False, read_callback=None):
        self.set_print_output(print_output)
        self.read_callback = read_callback
        self._stopped = False

    def stop_the_process(self):
        self._stopped = True
        self._kill()

    def _kill(self):
        if self.proc and self.proc.returncode is None:
            self.proc.kill()  # SIGKILL on Unix, TerminateProcess on Windows
        proc.terminate()  # SIGTERM on Unix, CTRL-BREAK-ish on Windows

        try:
            # 2) Give it time to clean up
            run_coro_sync(asyncio.wait_for(proc.wait(), 0.1))
            return
        except asyncio.TimeoutError:
            pass
        proc.kill()  # SIGKILL / TerminateProcess

    async def run_subprocess(self, read_args=[]):

        if self.directory:
            dr = os.getcwd()
            os.chdir(self.directory)
        else:
            dr = None
        if self._stopped:
            raise asyncio.CancelledError()

        self.proc = await asyncio.create_subprocess_exec(
            *self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **self.kwargs
        )
        if self._stopped:
            self._kill()

        if dr:
            os.chdir(dr)

        exception = None

        patch = self.outfile or self.print_output or self.read_callback

        def process(data, kind):
            utf8 = None
            if self.print_output:
                utf8 = data.decode("utf8")
                print(utf8, end="")
            if self.outfile:
                self.outfile.write(data)
            if self.read_callback:
                if utf8 is None:
                    utf8 = data.decode("utf8")
                self.read_callback(utf8, kind)

        def replace_feed_data(stream_reader, kind):
            nonlocal exception
            if stream_reader._buffer:
                process(stread_reader._buffer, kind)

            def feed_data(data):
                nonlocal exception
                try:
                    process(data, kind)
                except Exception as e:
                    if not exception:
                        exception = e
                    raise
                fd(data)

            fd = stream_reader.feed_data
            stream_reader.feed_data = feed_data

        if self.outfile:
            replace_feed_data(self.proc.stdout, "out")
            replace_feed_data(self.proc.stderr, "err")

        try:
            out = await asyncio.gather(
                self.read_output(self.proc.stdout, *read_args),
                self.read_error(self.proc.stderr, *read_args),
                self.proc.wait(),
            )
        finally:
            if self.outfile:
                self.outfile.close()

        if exception:
            raise exception

        return self.result(*out)

    def result(self, output, error, wait):
        """This function is for postprocessing the results.

        It is intended to be predefined in the descendants

        Parameters
        ----------

        output: mixed
          Result of the self.read_output

        error: mixed
          Result of the self.read_error

        wait: int
          The process return value

        Return
        ------

        out: mixed
          Currently, the tuple (output, error, return_code) is returned, however,
          subclasses can return anything they want.
        """

        if wait != 0:
            e = ValueError(f"The process ended with return value {wait}")
            e.output = output
            e.error = error
            raise e
        return output, error, wait

    def run_async(self, cmd, outfile, directory=None, read_args=[], **kwargs):
        self.cmd = cmd
        self.outfile = outfile
        self.directory = directory
        self.kwargs = kwargs
        return self.run_subprocess(read_args)

    def _run(self, cmd, outfile, directory=None, read_args=[], **kwargs):
        coro = self.run_async(cmd, outfile, directory, read_args, **kwargs)
        return run_coro_sync(coro)

    async def read_error(self, stderr, *args):
        while True:
            line = await stderr.readline()
            if not line:
                return
            print(line.decode("utf8"))

    async def read_output(self, stdout, *args):
        raise NotImplementedError(
            "Please, redefine ProcessOuputParser.read_output coroutine"
        )

    def set_print_output(self, print_output):
        self.print_output = print_output

    @maybeclassmethod
    def read_from_file(self, cls, output, error=None, read_args=[], return_code=0):
        """
        Synchronous wrapper for reading output/error from files.
        """
        if not self:
            self = cls()

        return run_coro_sync(
            self.read_output_file(output, error, read_args, return_code)
        )

    async def read_output_file(self, output, error=None, read_args=[], return_code=0):
        """
        Async version: reads output and error from files using AsyncioFileReader.
        """

        async def read_output_file():
            with AsyncioFileReader(output) as air:
                return await self.read_output(air, *read_args)

        async def read_error_file():
            if not error:
                return None
            with AsyncioFileReader(error) as air:
                return await self.read_error(air, *read_args)

        out_result, err_result = await asyncio.gather(
            read_output_file(), read_error_file()
        )

        return self.result(out_result, err_result, return_code)


class AsyncioFileReader:
    """File reader that mimics asyncio StreamReader.
    In fact, it is synchronous, but it offers the
    same interface as the asyncio StreamReader.
    """

    def __init__(self, filename, buffersize=8192):
        self.file = None  # just to set it in case the open fail
        self.file = open(filename, "rb") if isinstance(filename, str) else filename
        self.buffersize = buffersize

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

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

        # https://stackoverflow.com/questions/14128763/how-to-find-the-overlap-between-2-sequences-and-return-it
        shifts = self.separator_shifts(sep)

        # do the actual search
        # read the bytes into the array to avoid incrementing the array one by one
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


def run_coro_sync(coro_func):
    """
    Run an async coroutine from a synchronous context,
    even if an event loop is already running.
    `coro_func` must be an async function, not a coroutine object.
    """

    # Is there a running event loop?
    try:
        asyncio.get_running_loop()
        in_thread = True
    except RuntimeError:
        in_thread = False

    if not in_thread:
        return asyncio.run(coro_func)

    # YES → must run in another thread
    result = None
    exc = None

    def runner():
        nonlocal result, exc
        try:
            result = asyncio.run(coro_func)
        except Exception as e:
            exc = e

    t = Thread(target=runner)
    t.start()
    t.join()

    if exc:
        raise exc
    return result
