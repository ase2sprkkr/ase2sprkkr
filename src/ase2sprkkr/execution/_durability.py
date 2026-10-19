"""Crash-durable filesystem primitives for persistent execution state.

These helpers deliberately do not hide failed POSIX durability barriers.  A
caller that commits lifecycle state must either complete every requested
``fsync`` or receive an exception instead of being told that the transition is
durable.  Python does not expose an equivalent directory flush on Windows, so
there file contents are synchronized but directory-entry crash durability is a
weaker platform contract.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from uuid import uuid4


def fsync_directory(directory):
    """Synchronize one directory namespace where the platform supports it.

    POSIX errors deliberately propagate.  Windows has no portable directory
    ``fsync`` through Python, so callers still get file-data synchronization and
    atomic replacement, but not the stronger directory-entry persistence
    guarantee available on POSIX filesystems.
    """
    if os.name == "nt":
        return False

    directory = Path(directory)
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(directory, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return True


def durable_mkdir(path, *, parents=False, exist_ok=False, mode=0o777):
    """Create a directory and durably publish its path where supported.

    The final parent is synchronized even if the directory already exists.  It
    may have been created immediately before this call by infrastructure such
    as a lock helper, so existence alone is not proof that the directory entry
    crossed a persistence barrier.
    """
    path = Path(path)
    missing = []
    if parents:
        current = path
        while not current.exists():
            missing.append(current)
            parent = current.parent
            if parent == current:
                break
            current = parent
    elif not path.exists():
        missing.append(path)

    path.mkdir(parents=parents, exist_ok=exist_ok, mode=mode)
    if missing:
        for created in reversed(missing):
            fsync_directory(created.parent)
    else:
        fsync_directory(path.parent)
    return path


def durable_replace(source, destination):
    """Atomically replace *destination* and durably publish the rename."""
    source = Path(source)
    destination = Path(destination)
    os.replace(source, destination)
    fsync_directory(destination.parent)
    if source.parent != destination.parent:
        fsync_directory(source.parent)


def durable_write_text(path, text, *, encoding="utf-8", mode=None):
    """Atomically replace text after synchronizing contents and its new name."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    descriptor = None
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        descriptor = os.open(temporary, flags, 0o666 if mode is None else mode)
        if mode is not None:
            # Preserve the old best-effort permission contract.  Supplying the
            # mode to os.open prevents a wider-permission window before chmod,
            # while filesystems that reject chmod must not turn an otherwise
            # durable state write into a lifecycle failure.
            try:
                os.chmod(temporary, mode)
            except OSError:
                pass
        with os.fdopen(descriptor, "w", encoding=encoding) as stream:
            descriptor = None
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        durable_replace(temporary, path)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def durable_write_json(path, data, *, mode=None, **json_options):
    """Serialize JSON and commit it with :func:`durable_write_text`."""
    durable_write_text(
        path,
        json.dumps(data, **json_options),
        mode=mode,
    )


def sync_file(path):
    """Synchronize one existing file and the directory entry naming it."""
    path = Path(path)
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    fsync_directory(path.parent)


def sync_stream(stream):
    """Synchronize data written through an already open stream."""
    stream.flush()
    os.fsync(stream.fileno())
