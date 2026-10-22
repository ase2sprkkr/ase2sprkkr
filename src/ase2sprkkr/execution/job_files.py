"""File operations associated with one run's authoritative FileManifest."""

from __future__ import annotations

import os
from dataclasses import replace
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from uuid import uuid4

from ..common.file_versions import FileState, FileVersion, file_digest
from ._durability import sync_file
from .base import RunState
from .transport import expand_local_candidates, has_file_magic


class HistoricalFileMismatchError(OSError):
    """Restored contents differ from a run's captured SHA-256 digest."""


class ArtifactAvailabilityState(str, Enum):
    """Whether this run's declared artifacts are locally usable."""

    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"
    SKIPPED = "skipped"


_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _unique(values):
    return tuple(dict.fromkeys(value for value in values if value is not None))


class JobFiles:
    """File handling for a single persisted job, without parallel file state."""

    def __init__(self, record):
        self.record = record

    @property
    def manifest(self):
        return self.record.files

    @staticmethod
    def _concrete_result_files(manifests):
        paths = []
        seen = set()
        for manifest in manifests:
            for path in expand_local_candidates(manifest):
                path = os.path.abspath(path)
                if path in seen or not os.path.isfile(path):
                    continue
                seen.add(path)
                paths.append(path)
        return tuple(paths)

    @classmethod
    def _result_files_snapshot(cls, manifests, *, sync=False):
        versions = []
        for path in cls._concrete_result_files(manifests):
            if sync:
                sync_file(path)
            versions.append(FileVersion.capture(path, require_digest=True))
        return tuple(sorted(versions, key=lambda value: os.fspath(value.path)))

    @classmethod
    def _sync_result_files(cls, manifests):
        for path in cls._concrete_result_files(manifests):
            sync_file(path)

    @classmethod
    def _file_is_current(cls, result_files, path, versions):
        if result_files is None or path is None:
            return False
        path = os.path.abspath(path)
        if result_files.captured_at is None:
            return os.path.isfile(path)
        for version in versions:
            if os.path.abspath(os.fspath(version.path)) == path:
                return version.is_current()
        return False

    @staticmethod
    def _missing_result_files(result_files):
        if result_files is None:
            return ()
        if result_files.captured_at is not None:
            return tuple(
                os.fspath(item.path)
                for item in result_files.all_versions
                if item.state() is FileState.MISSING
            )

        missing = []
        for manifest in result_files.manifests:
            if has_file_magic(manifest):
                if not expand_local_candidates(manifest):
                    missing.append(manifest)
            elif not os.path.isfile(manifest):
                missing.append(manifest)
        return _unique(missing)

    @classmethod
    def _capture_result_files(cls, result_files):
        if result_files is None or result_files.captured_at is not None:
            return result_files
        from .file_availability import discover_local_files
        captured = replace(
            result_files,
            input_versions=cls._result_files_snapshot(result_files.input_files),
            output_versions=cls._result_files_snapshot(
                result_files.output_files, sync=True
            ),
            captured_at=_now(),
        )
        return discover_local_files(captured)

    def input_file_is_current(self):
        result_files = self.manifest
        return self._file_is_current(
            result_files,
            None if result_files is None else result_files.input_file,
            () if result_files is None else result_files.input_versions,
        )

    def output_file_is_current(self):
        result_files = self.manifest
        return self._file_is_current(
            result_files,
            None if result_files is None else result_files.output_file,
            () if result_files is None else result_files.output_versions,
        )

    def update_artifact_availability(self, store, path, state, message=None):
        state = ArtifactAvailabilityState(state)

        def update(latest):
            return replace(latest, updated_at=_now(), files=latest.files.with_file_availability(
                path, state.value, message,
            ))

        return store.update(self.record, update)

    def retrieve_result_files(self, store, executor=None):
        record = self.record
        if record.state not in _FINAL_STATES:
            raise RuntimeError(f"Run {record.id} is not finished")

        missing = self._missing_result_files(record.files)
        if not missing or record.execution.backend_id is None:
            return record

        if executor is None or not executor.supports_result_file_retrieval:
            return record

        versions = {
            os.path.abspath(os.fspath(version.path)): version
            for version in (() if record.files is None else record.files.all_versions)
        }
        for path in missing:
            path = os.path.abspath(os.fspath(path))
            version = versions.get(path)
            staged = (
                Path(path).with_name(f".{Path(path).name}.{uuid4().hex}.tmp")
                if version is not None and version.digest is not None else None
            )
            try:
                if staged is None:
                    executor.fetch_result_files(record, (path,))
                else:
                    executor.fetch_result_files(
                        record, (path,), destinations={path: os.fspath(staged)},
                    )
                    if not staged.is_file():
                        raise FileNotFoundError(
                            f"Historical artifact was not retrieved: {path}"
                        )
                    if file_digest(staged) != version.digest:
                        raise HistoricalFileMismatchError(
                            f"Historical result files have been overwritten: {path}"
                        )
                    os.replace(staged, path)
            except OSError as exc:
                self.update_artifact_availability(
                    store, path, ArtifactAvailabilityState.FAILED, str(exc),
                )
                raise
            finally:
                if staged is not None:
                    staged.unlink(missing_ok=True)
            refreshed = store.get(record.id)
            if path in self._missing_result_files(refreshed.files):
                record = JobFiles(refreshed).update_artifact_availability(
                    store, path, ArtifactAvailabilityState.FAILED,
                    "Historical artifact is still missing",
                )
            else:
                record = self.update_artifact_availability(
                    store,                     refreshed, path, ArtifactAvailabilityState.AVAILABLE,
                )
        return record


