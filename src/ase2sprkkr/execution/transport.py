"""Low-level asynchronous command and file transports."""

from __future__ import annotations

import asyncio
import os
import posixpath
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Protocol


class RunningCommand(Protocol):
    """A running command with independently readable output streams."""

    stdout: asyncio.StreamReader
    stderr: asyncio.StreamReader
    returncode: Optional[int]

    async def wait(self) -> int: ...

    def terminate(self) -> None: ...


class Transport:
    """Run commands and move files on one target machine.

    Transport operations are asynchronous because remote execution is mostly
    I/O. ``run()`` waits for a command and captures its output, while
    ``start()`` returns a process-like object with native asyncio stdout and
    stderr streams.
    """

    async def run(self, argv, *, cwd=None, input=None, check=True):
        raise NotImplementedError

    async def start(self, argv, *, cwd=None, stdin=None) -> RunningCommand:
        raise NotImplementedError

    async def put_file(self, source, destination):
        raise NotImplementedError

    async def get_file(self, source, destination):
        raise NotImplementedError

    async def read_file(self, path):
        raise NotImplementedError


async def _completed_process(process, argv, input=None, check=True):
    stdout, stderr = await process.communicate(input)
    result = subprocess.CompletedProcess(tuple(argv), process.returncode, stdout, stderr)
    if check and result.returncode:
        raise subprocess.CalledProcessError(
            result.returncode,
            result.args,
            output=result.stdout,
            stderr=result.stderr,
        )
    return result


class LocalTransport(Transport):
    """Transport backed by asyncio subprocesses and the local filesystem."""

    async def start(self, argv, *, cwd=None, stdin=None):
        return await asyncio.create_subprocess_exec(
            *(os.fspath(arg) for arg in argv),
            cwd=cwd,
            stdin=stdin,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )

    async def run(self, argv, *, cwd=None, input=None, check=True):
        argv = tuple(os.fspath(arg) for arg in argv)
        process = await asyncio.create_subprocess_exec(
            *argv,
            cwd=cwd,
            stdin=asyncio.subprocess.PIPE if input is not None else None,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        return await _completed_process(process, argv, input=input, check=check)

    async def put_file(self, source, destination):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(None, shutil.copyfile, source, destination)

    async def get_file(self, source, destination):
        await self.put_file(source, destination)

    async def read_file(self, path):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, Path(path).read_bytes)


class SshTransport(Transport):
    """Transport using the system ``ssh`` client."""

    def __init__(self, host, *, ssh="ssh", options=()):
        self.host = host
        self.ssh = ssh
        self.options = tuple(options)

    def _command(self, argv, cwd=None):
        command = shlex.join([os.fspath(arg) for arg in argv])
        if cwd is not None:
            command = "cd {} && exec {}".format(
                shlex.quote(os.fspath(cwd)), command
            )
        return [self.ssh, *self.options, self.host, command]

    async def start(self, argv, *, cwd=None, stdin=None):
        return await asyncio.create_subprocess_exec(
            *self._command(argv, cwd),
            stdin=stdin,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )

    async def run(self, argv, *, cwd=None, input=None, check=True):
        command = self._command(argv, cwd)
        process = await asyncio.create_subprocess_exec(
            *command,
            stdin=asyncio.subprocess.PIPE if input is not None else None,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        return await _completed_process(process, command, input=input, check=check)

    async def put_file(self, source, destination):
        destination = os.fspath(destination)
        parent = posixpath.dirname(destination)
        if parent:
            await self.run(("mkdir", "-p", "--", parent))

        command = self._command(("sh", "-c", 'cat > "$1"', "sh", destination))
        with open(source, "rb") as source_file:
            process = await asyncio.create_subprocess_exec(
                *command,
                stdin=source_file,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            await _completed_process(process, command)

    async def get_file(self, source, destination):
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        command = self._command(("cat", "--", os.fspath(source)))
        with destination.open("wb") as destination_file:
            process = await asyncio.create_subprocess_exec(
                *command,
                stdout=destination_file,
                stderr=asyncio.subprocess.PIPE,
            )
            _, stderr = await process.communicate()
            if process.returncode:
                raise subprocess.CalledProcessError(
                    process.returncode, command, stderr=stderr
                )

    async def read_file(self, path):
        return (await self.run(("cat", "--", os.fspath(path)))).stdout
