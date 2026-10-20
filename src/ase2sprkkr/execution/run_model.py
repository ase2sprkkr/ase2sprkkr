"""Persistent, non-redundant identity and file layout of one execution.

Submission intent and lifecycle are properties of a RunRecord, not backend
references.  This module contains only serializable value objects.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
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
        submitted = replace(self, backend_id=backend_id)
        if self.backend_id is not None and self.backend_id != submitted.backend_id:
            raise ValueError("Cannot replace an assigned backend identity")
        return submitted



class ArtifactAvailabilityState(str, Enum):
    """Availability of a declared local artifact."""

    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass(frozen=True)
class FileEntry:
    """Location of one file, with no assumptions about its content."""

    local_path: Optional[str] = None
    remote_path: Optional[str] = None
    availability_state: Optional[str] = None
    availability_message: Optional[str] = None
    source_pattern: Optional[str] = None

    def __post_init__(self):
        import os
        if self.source_pattern is not None:
            object.__setattr__(self, "source_pattern", os.fspath(self.source_pattern))

        for name in ("local_path", "remote_path"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, os.fspath(value))
        if self.local_path is None and self.remote_path is None:
            raise ValueError("A file entry must have at least one location")
        if self.availability_state is not None:
            try:
                availability = ArtifactAvailabilityState(self.availability_state)
            except ValueError as exc:
                raise ValueError("Invalid per-file availability state") from exc
            # Persist JSON-compatible strings, not implementation-specific enums.
            object.__setattr__(self, "availability_state", availability.value)
        if (
            self.availability_message is not None
            and self.availability_state != ArtifactAvailabilityState.FAILED
        ):
            raise ValueError("An availability error belongs to a failed file")


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


class FileSyncPolicy(str, Enum):
    """Policy controlling transfers between local and remote workspaces."""

    CHANGED = "changed"
    ALWAYS = "always"
    MISSING = "missing"
    NEVER = "never"

    def should_transfer(self, source_digest, destination_digest):
        """Apply the direction-independent transfer rule to file digests.

        The caller handles missing mandatory files separately; a missing
        source is never transferable, including with ALWAYS.
        """
        if self == FileSyncPolicy.NEVER or source_digest is None:
            return False
        if self == FileSyncPolicy.ALWAYS:
            return True
        if self == FileSyncPolicy.MISSING:
            return destination_digest is None
        return source_digest != destination_digest


@dataclass(frozen=True)
class FileMapping:
    """Resolved mapping and transfer policy of this run's file manifest."""

    local_root: str
    remote_root: str
    input_sync: str = FileSyncPolicy.CHANGED.value
    output_sync: str = FileSyncPolicy.CHANGED.value

    def __post_init__(self):
        import os

        object.__setattr__(self, "local_root", os.fspath(self.local_root))
        object.__setattr__(self, "remote_root", os.fspath(self.remote_root))
        try:
            input_sync = FileSyncPolicy(self.input_sync)
            output_sync = FileSyncPolicy(self.output_sync)
        except ValueError as exc:
            raise ValueError("Invalid file synchronization policy") from exc
        object.__setattr__(self, "input_sync", input_sync.value)
        object.__setattr__(self, "output_sync", output_sync.value)


@dataclass(frozen=True)
class FileManifest:
    """All paths and cache baselines required for transfer and reconnect."""

    inputs: tuple[FileEntry, ...] = ()
    outputs: tuple[FileEntry, ...] = ()
    caches: tuple[CacheEntry, ...] = ()
    stdout: Optional[FileEntry] = None
    stderr: Optional[FileEntry] = None
    mapping: Optional[FileMapping] = None
    input_versions: tuple[FileVersion, ...] = ()
    output_versions: tuple[FileVersion, ...] = ()
    captured_at: Optional[str] = None

    @classmethod
    def from_spec(cls, spec):
        """Build the local file manifest for an execution specification."""
        import os

        from .base import file_path

        directory = os.path.abspath(os.fspath(spec.directory))

        def paths(*values):
            return tuple(dict.fromkeys(
                path for value in values
                if (path := file_path(value, directory)) is not None
            ))

        inputs = paths(getattr(spec, "input_file", None), *getattr(spec, "input_files", ()))
        stdout = file_path(getattr(spec, "output_file", None), directory)
        stderr = file_path(getattr(spec, "error_file", None), directory)
        streams = {stdout, stderr}
        outputs = (
            path for path in paths(*getattr(spec, "output_files", ())) if path not in streams
        )
        caches = paths(*getattr(spec, "cache_files", ()))
        return cls(
            inputs=tuple(FileEntry(local_path=x) for x in inputs),
            outputs=tuple(FileEntry(local_path=x) for x in outputs),
            caches=tuple(CacheEntry(FileEntry(local_path=x)) for x in caches),
            stdout=None if stdout is None else FileEntry(local_path=stdout),
            stderr=None if stderr is None else FileEntry(local_path=stderr),
        ).initialize_availability()

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
        if isinstance(self.mapping, dict):
            object.__setattr__(self, "mapping", FileMapping(**self.mapping))
        for name in ("input_versions", "output_versions"):
            object.__setattr__(
                self, name,
                tuple(
                    version if isinstance(version, FileVersion) else FileVersion(**version)
                    for version in getattr(self, name)
                ),
            )
    @property
    def artifacts(self):
        """Output streams, declared outputs and caches (one status per local file)."""
        entries = (self.stdout, self.stderr, *self.outputs, *(x.file for x in self.caches))
        unique = {}
        for entry in entries:
            if entry is not None and entry.local_path is not None:
                unique[entry.local_path] = entry
        return tuple(unique.values())

    @property
    def failed_artifacts(self):
        return tuple(x for x in self.artifacts if x.availability_state == "failed")

    @property
    def availability_state(self):
        """Derived summary for callers; never stored independently."""
        statuses = tuple(x.availability_state for x in self.artifacts)
        if "failed" in statuses:
            return "failed"
        if not statuses or None in statuses or "pending" in statuses:
            return "pending"
        if all(x == "skipped" for x in statuses):
            return "skipped"
        return "available"

    @property
    def availability_message(self):
        failures = self.failed_artifacts
        return None if not failures else failures[0].availability_message

    def with_remote_stream_paths(self, *, stdout=None, stderr=None):
        """Update recovered stream locations without changing artifact identity.

        Missing or empty recovered paths do not replace persisted locations.
        A stream that was not recorded remains absent.
        """
        changes = {}
        for name, remote in (("stdout", stdout), ("stderr", stderr)):
            current = getattr(self, name)
            if remote and current is not None:
                changes[name] = replace(current, remote_path=remote)
        return replace(self, **changes) if changes else self

    def initialize_availability(self):
        """Mark declared downloadable files pending without marking inputs."""
        return self._map_artifacts(
            lambda item: replace(item, availability_state="pending")
            if item.availability_state is None else item
        )

    def _map_artifacts(self, transform):
        return replace(
            self,
            stdout=None if self.stdout is None else transform(self.stdout),
            stderr=None if self.stderr is None else transform(self.stderr),
            outputs=tuple(transform(x) for x in self.outputs),
            caches=tuple(replace(x, file=transform(x.file)) for x in self.caches),
        )

    def with_file_availability(
        self, local_path, state, message=None, *,
        remote_path=None, kind="output", source_pattern=None,
    ):
        """Update one artifact, adding each wildcard match as it is discovered.

        The pattern declaration retains its own independent status. A child's
        source_pattern identifies the declaration, including after restart.
        Duplicate stream/output declarations of one path are updated together.
        """
        import os
        path = os.fspath(local_path)
        matches = any(x.local_path == path for x in self.artifacts)
        if matches:
            return self._map_artifacts(
                lambda item: replace(
                    item,
                    availability_state=state,
                    availability_message=message,
                ) if item.local_path == path else item
            )
        file = FileEntry(
            local_path=path, remote_path=remote_path,
            availability_state=state, availability_message=message,
            source_pattern=source_pattern,
        )
        if kind == "output":
            return replace(self, outputs=(*self.outputs, file))
        if kind == "cache":
            return replace(self, caches=(*self.caches, CacheEntry(file)))
        raise ValueError(f"Unknown artifact kind: {kind!r}")

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
        streams = (self.stdout, self.stderr)
        entries = tuple(entry for entry in streams if entry is not None)
        return tuple(dict.fromkeys(
            entry.local_path
            for entry in (*entries, *self.outputs)
            if entry.local_path is not None
        ))

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
