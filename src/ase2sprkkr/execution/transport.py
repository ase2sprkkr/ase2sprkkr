"""Low-level asynchronous command and file transports."""

from __future__ import annotations

import asyncio
import fnmatch
import glob
import os
import posixpath
import shlex
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path, PurePosixPath

from ..common.file_versions import file_digest
from typing import Optional, Protocol
from uuid import uuid4


class RunningCommand(Protocol):
    """A running command with independently readable output streams."""

    stdout: asyncio.StreamReader
    stderr: asyncio.StreamReader
    returncode: Optional[int]

    async def wait(self) -> int: ...

    def terminate(self) -> None: ...


def _sources(sources):
    if isinstance(sources, (str, bytes, os.PathLike)):
        return (sources,)
    return tuple(sources)


def has_file_magic(path):
    """Return whether *path* contains shell-style wildcard characters."""
    return glob.has_magic(os.fspath(path))


def expand_local_candidates(path):
    """Enumerate path candidates; preserve exact names even when absent.

    No root validation or symlink resolution happens here. The operation
    consuming the paths is responsible for enforcing its own safety policy.
    """
    path = os.fspath(path)
    if has_file_magic(path):
        return tuple(sorted(glob.glob(path, recursive=True)))
    return (path,)


def _local_pattern(source, root):
    root = Path(root).expanduser().resolve()
    source = os.fspath(source)
    path = Path(source).expanduser()
    if not path.is_absolute():
        path = root / path
    pattern = os.path.abspath(os.fspath(path))
    try:
        os.path.relpath(pattern, root)
    except ValueError as exc:
        raise ValueError(f"File pattern {source!r} is outside source root {root}") from exc
    relative = os.path.relpath(pattern, root)
    if relative == os.pardir or relative.startswith(os.pardir + os.sep):
        raise ValueError(f"File pattern {source!r} is outside source root {root}")
    return root, pattern


def expand_local_files(sources, root):
    """Expand local source paths/patterns below *root* into concrete files."""
    root = Path(root).expanduser().resolve()
    files = []
    seen = set()
    for source in _sources(sources):
        _, pattern = _local_pattern(source, root)
        for match in expand_local_candidates(pattern):
            path = Path(match)
            if not path.is_file():
                continue
            path = path.resolve()
            try:
                path.relative_to(root)
            except ValueError as exc:
                raise ValueError(
                    f"Expanded file {path} is outside source root {root}"
                ) from exc
            value = os.fspath(path)
            if value not in seen:
                seen.add(value)
                files.append(value)
    return tuple(files)


def _posix_relative_pattern(source, root):
    source = os.fspath(source)
    root = posixpath.normpath(os.fspath(root))
    if posixpath.isabs(source):
        pattern = posixpath.normpath(source)
        try:
            common = posixpath.commonpath((root, pattern))
        except ValueError as exc:
            raise ValueError(
                f"File pattern {source!r} is outside source root {root}"
            ) from exc
        if common != root:
            raise ValueError(f"File pattern {source!r} is outside source root {root}")
        relative = posixpath.relpath(pattern, root)
    else:
        relative = posixpath.normpath(source)
    if relative == ".." or relative.startswith("../"):
        raise ValueError(f"File pattern {source!r} is outside source root {root}")
    return relative


def _match_component(name, pattern):
    if name.startswith(".") and not pattern.startswith("."):
        return False
    return fnmatch.fnmatchcase(name, pattern)


def _match_relative_path(path, pattern):
    parts = PurePosixPath(path).parts
    patterns = PurePosixPath(pattern).parts

    @lru_cache(maxsize=None)
    def match(path_index, pattern_index):
        if pattern_index == len(patterns):
            return path_index == len(parts)
        current = patterns[pattern_index]
        if current == "**":
            # Match glob.glob(recursive=True): an implicit wildcard never
            # descends into a hidden directory or consumes a hidden file.
            return match(path_index, pattern_index + 1) or (
                path_index < len(parts)
                and not parts[path_index].startswith(".")
                and match(path_index + 1, pattern_index)
            )
        return (
            path_index < len(parts)
            and _match_component(parts[path_index], current)
            and match(path_index + 1, pattern_index + 1)
        )

    return match(0, 0)


def _literal_prefix(pattern):
    parts = PurePosixPath(pattern).parts
    prefix = []
    for part in parts:
        if has_file_magic(part) or part == "**":
            break
        prefix.append(part)
    return prefix


def _local_temporary_path(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination.with_name(f".{destination.name}.{uuid4().hex}.tmp")


class Transport:
    """Run commands and move files on one target machine.

    Transport operations are asynchronous because remote execution is mostly
    I/O. ``put_file()`` and ``get_file()`` are exact one-file primitives.
    ``put_files()`` and ``get_files()`` accept one path/pattern or an iterable
    and preserve each matched file's path relative to the supplied source root.
    Wildcards are always expanded on the source side of the transfer.
    """

    async def run(self, argv, *, cwd=None, input=None, check=True):
        raise NotImplementedError

    async def start(self, argv, *, cwd=None, stdin=None) -> RunningCommand:
        raise NotImplementedError

    async def put_file(self, source, destination):
        raise NotImplementedError

    async def get_file(self, source, destination):
        raise NotImplementedError

    async def put_files(self, sources, source_root, destination_root):
        raise NotImplementedError

    async def get_files(self, sources, source_root, destination_root):
        raise NotImplementedError

    async def expand_files(self, sources, root):
        """Expand paths/patterns on the target filesystem into concrete files."""
        raise NotImplementedError

    async def read_file(self, path):
        raise NotImplementedError

    async def file_digest(self, path):
        """Return the SHA-256 digest of *path*, or ``None`` if it does not exist."""
        raise NotImplementedError

    async def file_signature(self, path):
        """Return ``(size, mtime_ns)`` for *path*, or ``None`` if absent."""
        raise NotImplementedError



def local_file_signature(path):
    """Return ``(size, mtime_ns)`` for a local file, or ``None`` if absent."""
    try:
        stat = os.stat(path)
    except FileNotFoundError:
        return None
    return stat.st_size, stat.st_mtime_ns


def _parse_remote_file_signature(raw):
    """Parse GNU ``stat`` output into ``(size, mtime_ns)``."""
    size, seconds, timestamp = raw.decode("utf8").rstrip("\n").split("|", 2)
    fraction = ""
    if "." in timestamp:
        fraction = timestamp.split(".", 1)[1].split(None, 1)[0]
    nanoseconds = int((fraction + "000000000")[:9]) if fraction else 0
    return int(size), int(seconds) * 1_000_000_000 + nanoseconds


# Kept temporarily for code using the helper inside the execution package.


async def _completed_process(process, argv, input=None, check=True):
    try:
        stdout, stderr = await process.communicate(input)
    except BaseException:
        # Cancelling a command must not abandon its subprocess (or leave an
        # unreaped zombie). This also handles failures during pipe I/O.
        if process.returncode is None:
            try:
                process.kill()
            except ProcessLookupError:
                pass
        await process.wait()
        raise
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
        source = os.fspath(source)
        signature = local_file_signature(source)
        if signature is None:
            raise FileNotFoundError(source)
        destination = Path(destination)
        temporary = _local_temporary_path(destination)
        loop = asyncio.get_running_loop()
        try:
            await loop.run_in_executor(None, shutil.copyfile, source, temporary)
            if local_file_signature(source) != signature:
                raise RuntimeError(f"Local file changed while copying: {source}")
            stat = temporary.stat()
            os.utime(temporary, ns=(stat.st_atime_ns, signature[1]))
            os.replace(temporary, destination)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    async def get_file(self, source, destination):
        await self.put_file(source, destination)

    async def expand_files(self, sources, root):
        return expand_local_files(sources, root)

    async def _transfer_files(self, paths, source_root, destination_root, transfer):
        """Preserve relative layout when copying selected local files."""
        source_root = Path(source_root).expanduser().resolve()
        destination_root = Path(destination_root).expanduser().resolve()
        transfers = []
        for source in paths:
            destination = destination_root / Path(source).relative_to(source_root)
            await transfer(source, destination)
            transfers.append((source, os.fspath(destination)))
        return tuple(transfers)

    async def put_files(self, sources, source_root, destination_root):
        paths = expand_local_files(sources, source_root)
        return await self._transfer_files(paths, source_root, destination_root, self.put_file)

    async def get_files(self, sources, source_root, destination_root):
        paths = await self.expand_files(sources, source_root)
        return await self._transfer_files(paths, source_root, destination_root, self.get_file)

    async def read_file(self, path):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, Path(path).read_bytes)

    async def file_digest(self, path):
        loop = asyncio.get_running_loop()
        return await asyncio.to_thread(file_digest, path, missing_ok=True)

    async def file_signature(self, path):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, local_file_signature, path)


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
        source = os.fspath(source)
        signature = local_file_signature(source)
        if signature is None:
            raise FileNotFoundError(source)

        destination = os.fspath(destination)
        parent = posixpath.dirname(destination)
        if parent:
            await self.run(("mkdir", "-p", "--", parent))
        temporary = posixpath.join(
            parent or ".",
            f".{posixpath.basename(destination)}.{uuid4().hex}.tmp",
        )

        command = self._command(("sh", "-c", 'cat > "$1"', "sh", temporary))
        try:
            with open(source, "rb") as source_file:
                process = await asyncio.create_subprocess_exec(
                    *command,
                    stdin=source_file,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )
                await _completed_process(process, command)
            if local_file_signature(source) != signature:
                raise RuntimeError(f"Local file changed while uploading: {source}")
            seconds, nanoseconds = divmod(signature[1], 1_000_000_000)
            await self.run(
                (
                    "touch",
                    "-d",
                    f"@{seconds}.{nanoseconds:09d}",
                    "--",
                    temporary,
                )
            )
            await self.run(("mv", "-f", "--", temporary, destination))
        except BaseException:
            await self.run(("rm", "-f", "--", temporary), check=False)
            raise

    async def get_file(self, source, destination):
        source = os.fspath(source)
        signature = await self.file_signature(source)
        if signature is None:
            raise FileNotFoundError(source)

        destination = Path(destination)
        temporary = _local_temporary_path(destination)
        command = self._command(("cat", "--", source))
        try:
            with temporary.open("wb") as destination_file:
                process = await asyncio.create_subprocess_exec(
                    *command,
                    stdout=destination_file,
                    stderr=asyncio.subprocess.PIPE,
                )
                await _completed_process(process, command)
            if await self.file_signature(source) != signature:
                raise RuntimeError(f"Remote file changed while downloading: {source}")
            stat = temporary.stat()
            os.utime(temporary, ns=(stat.st_atime_ns, signature[1]))
            os.replace(temporary, destination)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    async def expand_files(self, sources, root):
        root = posixpath.normpath(os.fspath(root))
        files = []
        seen = set()
        for source in _sources(sources):
            pattern = _posix_relative_pattern(source, root)
            if not has_file_magic(pattern):
                path = posixpath.join(root, pattern)
                if await self.file_signature(path) is not None and path not in seen:
                    seen.add(path)
                    files.append(path)
                continue

            prefix = _literal_prefix(pattern)
            search_root = posixpath.join(root, *prefix) if prefix else root
            exists = await self.run(
                ("sh", "-c", '[ -d "$1" ] || exit 3', "sh", search_root),
                check=False,
            )
            if exists.returncode == 3:
                continue
            if exists.returncode:
                raise subprocess.CalledProcessError(
                    exists.returncode, exists.args,
                    output=exists.stdout, stderr=exists.stderr,
                )
            result = await self.run(("find", search_root, "-type", "f", "-print0"))
            for item in result.stdout.decode("utf8").split("\0"):
                if not item:
                    continue
                relative = posixpath.relpath(item, root)
                if _match_relative_path(relative, pattern) and item not in seen:
                    seen.add(item)
                    files.append(item)
        return tuple(sorted(files))

    @staticmethod
    async def _pipe(source, destination):
        try:
            while True:
                chunk = await source.read(1024 * 1024)
                if not chunk:
                    break
                destination.write(chunk)
                await destination.drain()
        finally:
            destination.close()

    async def _transfer_tar(self, producer_command, consumer_command):
        """Stream an archive between two subprocesses and check both results."""
        consumer = await asyncio.create_subprocess_exec(
            *consumer_command,
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.DEVNULL,
            stderr=asyncio.subprocess.PIPE,
        )
        producer = None
        try:
            producer = await asyncio.create_subprocess_exec(
                *producer_command,
                stdin=asyncio.subprocess.DEVNULL,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            async def pipe():
                try:
                    await self._pipe(producer.stdout, consumer.stdin)
                except (BrokenPipeError, ConnectionResetError):
                    # If the extractor fails early, the archive producer can
                    # otherwise remain blocked on a full stdout pipe.
                    if producer.returncode is None:
                        producer.kill()

            _, producer_error, consumer_error, _, _ = await asyncio.gather(
                pipe(),
                producer.stderr.read(),
                consumer.stderr.read(),
                producer.wait(),
                consumer.wait(),
            )
            for process, command, stderr in (
                (producer, producer_command, producer_error),
                (consumer, consumer_command, consumer_error),
            ):
                if process.returncode:
                    raise subprocess.CalledProcessError(
                        process.returncode, command, stderr=stderr,
                    )
        finally:
            for process in (producer, consumer):
                if process is not None and process.returncode is None:
                    process.kill()
            await asyncio.gather(
                *(process.wait() for process in (producer, consumer)
                  if process is not None),
            )

    async def put_files(self, sources, source_root, destination_root):
        source_root = Path(source_root).expanduser().resolve()
        destination_root = posixpath.normpath(os.fspath(destination_root))
        files = expand_local_files(sources, source_root)
        if not files:
            return ()

        relatives = [os.fspath(Path(source).relative_to(source_root)) for source in files]
        await self.run(("mkdir", "-p", "--", destination_root))
        producer = ("tar", "-C", os.fspath(source_root), "-cf", "-", "--", *relatives)
        consumer = self._command(("tar", "-C", destination_root, "-xf", "-"))
        await self._transfer_tar(producer, consumer)
        return tuple(
            (source, posixpath.join(destination_root, *Path(relative).parts))
            for source, relative in zip(files, relatives)
        )

    async def get_files(self, sources, source_root, destination_root):
        source_root = posixpath.normpath(os.fspath(source_root))
        destination_root = Path(destination_root).expanduser().resolve()
        files = await self.expand_files(sources, source_root)
        if not files:
            return ()

        relatives = [posixpath.relpath(source, source_root) for source in files]
        destination_root.mkdir(parents=True, exist_ok=True)
        producer = self._command(
            ("tar", "-C", source_root, "-cf", "-", "--", *relatives)
        )
        consumer = ("tar", "-C", os.fspath(destination_root), "-xf", "-")
        await self._transfer_tar(producer, consumer)
        return tuple(
            (source, os.fspath(destination_root.joinpath(*PurePosixPath(relative).parts)))
            for source, relative in zip(files, relatives)
        )

    async def read_file(self, path):
        return (await self.run(("cat", "--", os.fspath(path)))).stdout

    async def _file_probe(self, path, command):
        """Run a remote file query; absent files are not command failures."""
        result = await self.run(
            ("sh", "-c", f'[ -f "$1" ] || exit 3; {command}', "sh", os.fspath(path)),
            check=False,
        )
        if result.returncode == 3:
            return None
        if result.returncode:
            raise subprocess.CalledProcessError(
                result.returncode,
                result.args,
                output=result.stdout,
                stderr=result.stderr,
            )
        return result.stdout

    async def file_digest(self, path):
        output = await self._file_probe(path, 'sha256sum -- "$1"')
        return None if output is None else output.decode("ascii").split(None, 1)[0]

    async def file_signature(self, path):
        output = await self._file_probe(path, 'stat -c "%s|%Y|%y" -- "$1"')
        return None if output is None else _parse_remote_file_signature(output)
