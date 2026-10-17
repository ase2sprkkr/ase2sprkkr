"""Persistent, non-redundant identity and file layout of one execution.

Submission intent and lifecycle are properties of a RunRecord, not backend
references.  This module contains only serializable value objects.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

from ..common.file_versions import FileState, FileVersion

from .endpoint import ExecutorEndpoint


@dataclass(frozen=True)
class WorkingDirectory:
    """Local working directory and its optional target-side counterpart."""

    local: str
    remote: Optional[str] = None

    def __post_init__(self):
        import os

        object.__setattr__(self, "local", os.fspath(self.local))
        if self.remote is not None:
            object.__setattr__(self, "remote", os.fspath(self.remote))

    @property
    def target(self):
        return self.remote if self.remote is not None else self.local


@dataclass(frozen=True)
class Execution:
    """Backend identity plus the directory in which the job executes."""

    endpoint: ExecutorEndpoint
    directory: WorkingDirectory
    backend_id: Optional[str] = None

    def __post_init__(self):
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        if isinstance(self.directory, dict):
            object.__setattr__(self, "directory", WorkingDirectory(**self.directory))
        if self.backend_id is not None:
            identifier = str(self.backend_id)
            if not identifier:
                raise ValueError("backend_id cannot be empty")
            object.__setattr__(self, "backend_id", identifier)


    def submitted(self, backend_id: str):
        """Return the execution bound to its immutable backend identity."""
        backend_id = str(backend_id)
        if not backend_id:
            raise ValueError("backend_id cannot be empty")
        if self.backend_id is not None and self.backend_id != backend_id:
            raise ValueError("Cannot replace an assigned backend identity")
        return replace(self, backend_id=backend_id)



@dataclass(frozen=True)
class FileEntry:
    """Location of one file, with no assumptions about its content."""

    local_path: Optional[str] = None
    remote_path: Optional[str] = None

    def __post_init__(self):
        import os

        for name in ("local_path", "remote_path"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, os.fspath(value))
        if self.local_path is None and self.remote_path is None:
            raise ValueError("A file entry must have at least one location")


@dataclass(frozen=True)
class CacheEntry:
    """Cache identity and its pre-execution remote signature."""

    file: FileEntry
    before: Optional[tuple] = None

    def __post_init__(self):
        if isinstance(self.file, dict):
            object.__setattr__(self, "file", FileEntry(**self.file))
        if self.before is not None:
            object.__setattr__(self, "before", tuple(self.before))


@dataclass(frozen=True)
class FileManifest:
    """All paths and cache baselines required for transfer and reconnect."""

    inputs: tuple[FileEntry, ...] = ()
    outputs: tuple[FileEntry, ...] = ()
    caches: tuple[CacheEntry, ...] = ()
    stdout: Optional[FileEntry] = None
    stderr: Optional[FileEntry] = None
    input_versions: tuple[FileVersion, ...] = ()
    output_versions: tuple[FileVersion, ...] = ()
    captured_at: Optional[str] = None
    availability_state: Optional[str] = None
    availability_message: Optional[str] = None

    def __post_init__(self):
        for name, cls in (("inputs", FileEntry), ("outputs", FileEntry), ("caches", CacheEntry)):
            items = getattr(self, name)
            object.__setattr__(
                self, name,
                tuple(item if isinstance(item, cls) else cls(**item) for item in items),
            )
        for name in ("stdout", "stderr"):
            item = getattr(self, name)
            if isinstance(item, dict):
                object.__setattr__(self, name, FileEntry(**item))
        for name in ("input_versions", "output_versions"):
            object.__setattr__(
                self, name,
                tuple(
                    version if isinstance(version, FileVersion) else FileVersion(**version)
                    for version in getattr(self, name)
                ),
            )
        if self.availability_state is not None:
            allowed = {"pending", "available", "failed"}
            if self.availability_state not in allowed:
                raise ValueError("Invalid artifact availability state")

    @property
    def input_file(self):
        return self.inputs[0].local_path if self.inputs else None

    @property
    def output_file(self):
        return self.stdout.local_path if self.stdout is not None else None

    @property
    def input_files(self):
        return tuple(entry.local_path for entry in self.inputs if entry.local_path is not None)

    @property
    def output_files(self):
        entries = (self.stdout,) if self.stdout is not None else ()
        return tuple(entry.local_path for entry in (*entries, *self.outputs) if entry.local_path is not None)

    @property
    def manifests(self):
        return tuple(dict.fromkeys((*self.input_files, *self.output_files)))

    @property
    def all_versions(self):
        by_path = {
            str(item.path): item
            for item in (*self.input_versions, *self.output_versions)
        }
        return tuple(by_path[key] for key in sorted(by_path))

    @property
    def state(self):
        if self.captured_at is None:
            return FileState.CURRENT
        states = tuple(version.state() for version in self.output_versions)
        if all(value is FileState.CURRENT for value in states):
            return FileState.CURRENT
        if FileState.OVERWRITTEN in states:
            return FileState.OVERWRITTEN
        if FileState.UNAVAILABLE in states:
            return FileState.UNAVAILABLE
        return FileState.MISSING
