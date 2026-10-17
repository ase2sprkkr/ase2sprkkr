"""Persistent run history and lifecycle orchestration."""

from __future__ import annotations

import asyncio
import glob
import json
import os
import warnings
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from threading import RLock
from typing import Any, Optional
from uuid import uuid4

from platformdirs import user_state_dir

from ..common.file_locks import FileLock
from ..common.file_versions import FileState, FileVersion
from ._durability import durable_mkdir, durable_write_json, sync_file
from .base import (
    Executor,
    FailurePhase,
    RunReference,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .endpoint import ExecutorEndpoint, endpoint_for_executor
from .pool import default_executor_pool
from .postprocessing import PostprocessingRunner


_RUN_RECORD_SCHEMA_VERSION = 2
_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}
_STATE_RANK = {
    RunState.NEW: 0,
    RunState.SUBMITTING: 1,
    RunState.QUEUED: 2,
    RunState.RUNNING: 3,
    RunState.DOWNLOADING: 4,
    RunState.POSTPROCESSING: 5,
    RunState.FINISHED: 6,
    RunState.FAILED: 6,
    RunState.CANCELLED: 6,
}
_UNCHANGED = object()
_NO_ARTIFACTS_FAILURES = {
    FailurePhase.SUBMISSION,
    FailurePhase.LAUNCH,
    FailurePhase.ARTIFACTS,
    FailurePhase.DOWNLOAD,
}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _timestamp(value):
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed


def _result_file_path(value, directory):
    if value is None:
        return None
    value = getattr(value, "name", value)
    try:
        path = os.fspath(value)
    except TypeError:
        return None
    if isinstance(path, int):
        return None
    if not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


def _unique(values):
    return tuple(dict.fromkeys(value for value in values if value is not None))


def _failure_phase(value):
    if value is None or isinstance(value, FailurePhase):
        return value
    return FailurePhase(value)


class CorruptRunRecord(RuntimeError):
    """A persisted run record cannot be decoded or validated."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Corrupt run record {self.path}: {message}")


class ArtifactAvailabilityState(str, Enum):
    """Whether this run's declared artifacts are locally usable."""

    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"


@dataclass(frozen=True)
class ResultFilesRecord:
    """Declared result files and immutable versions captured at completion."""

    input_file: Optional[str] = None
    output_file: Optional[str] = None
    input_files: tuple[str, ...] = ()
    output_files: tuple[str, ...] = ()
    input_versions: tuple[FileVersion, ...] = ()
    output_versions: tuple[FileVersion, ...] = ()
    captured_at: Optional[str] = None
    availability_state: Optional[ArtifactAvailabilityState] = None
    availability_message: Optional[str] = None

    def __post_init__(self):
        for name in ("input_file", "output_file"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, os.fspath(value))
        for name in ("input_files", "output_files"):
            object.__setattr__(
                self,
                name,
                tuple(os.fspath(value) for value in getattr(self, name)),
            )
        for name in ("input_versions", "output_versions"):
            object.__setattr__(
                self,
                name,
                tuple(
                    value if isinstance(value, FileVersion) else FileVersion(**value)
                    for value in getattr(self, name)
                ),
            )
        if self.availability_state is not None and not isinstance(
            self.availability_state, ArtifactAvailabilityState
        ):
            object.__setattr__(
                self,
                "availability_state",
                ArtifactAvailabilityState(self.availability_state),
            )

    @property
    def manifests(self):
        return _unique((*self.input_files, *self.output_files))

    @property
    def all_versions(self):
        return tuple(
            sorted(
                {
                    os.fspath(item.path): item
                    for item in (*self.input_versions, *self.output_versions)
                }.values(),
                key=lambda item: os.fspath(item.path),
            )
        )

    @property
    def state(self):
        if self.captured_at is None:
            return FileState.CURRENT
        states = tuple(version.state() for version in self.output_versions)
        if all(state is FileState.CURRENT for state in states):
            return FileState.CURRENT
        if FileState.OVERWRITTEN in states:
            return FileState.OVERWRITTEN
        if FileState.UNAVAILABLE in states:
            return FileState.UNAVAILABLE
        return FileState.MISSING


@dataclass(frozen=True)
class SubmissionKey:
    """Persistent duplicate-submission identity and recovery context."""

    endpoint: ExecutorEndpoint
    directory: str
    task: str
    dataset: Optional[str] = None
    context: Any = None

    def __post_init__(self):
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        object.__setattr__(self, "directory", os.fspath(self.directory))
        object.__setattr__(self, "task", str(self.task).upper())
        if self.dataset is not None:
            object.__setattr__(self, "dataset", str(self.dataset))

    def conflicts_with(self, other):
        if self.endpoint != other.endpoint or self.directory != other.directory:
            return False
        if self.task == "SCF" or other.task == "SCF":
            return True
        return self.task == other.task and self.dataset == other.dataset


class DuplicateSubmissionError(RuntimeError):
    """A persistent run would duplicate an active calculation namespace."""

    def __init__(self, record, conflicts):
        self.record = record
        self.conflicts = tuple(conflicts)
        super().__init__(
            "Run conflicts with {} active persistent run{}".format(
                len(self.conflicts), "" if len(self.conflicts) == 1 else "s"
            )
        )


@dataclass(frozen=True)
class RunRecord:
    """Persistent identity and authoritative lifecycle state of one run."""

    id: str
    state: RunState
    created_at: str
    updated_at: str
    task: str
    directory: str
    finished_at: Optional[str] = None
    dataset: Optional[str] = None
    endpoint: Optional[ExecutorEndpoint] = None
    submission: Optional[SubmissionKey] = None
    reference: Optional[RunReference] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    result_files: Optional[ResultFilesRecord] = None
    input_definition: Optional[str] = None
    failure_phase: Optional[FailurePhase] = None
    schema_version: int = _RUN_RECORD_SCHEMA_VERSION

    def __post_init__(self):
        object.__setattr__(self, "task", str(self.task).upper())
        if self.input_definition is not None:
            object.__setattr__(
                self, "input_definition", str(self.input_definition).upper()
            )
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        if isinstance(self.submission, dict):
            object.__setattr__(self, "submission", SubmissionKey(**self.submission))
        if isinstance(self.reference, dict):
            object.__setattr__(self, "reference", RunReference(**self.reference))
        if isinstance(self.result_files, dict):
            object.__setattr__(
                self, "result_files", ResultFilesRecord(**self.result_files)
            )
        if not isinstance(self.state, RunState):
            object.__setattr__(self, "state", RunState(self.state))
        if self.failure_phase is not None and not isinstance(
            self.failure_phase, FailurePhase
        ):
            object.__setattr__(
                self, "failure_phase", FailurePhase(self.failure_phase)
            )
        if self.schema_version != _RUN_RECORD_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported schema version {self.schema_version}; "
                f"expected {_RUN_RECORD_SCHEMA_VERSION}"
            )


class RunStore:
    """Persistent run history stored as date-sharded durable JSON files."""

    def __init__(self, directory=None):
        if directory is None:
            directory = Path(user_state_dir("ase2sprkkr")) / "runs"
        self.directory = Path(directory).expanduser()
        self._lock = RLock()
        self._warned_corrupt = set()

    @staticmethod
    def _bucket(created_at):
        try:
            return datetime.fromisoformat(created_at).date().isoformat()
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid run creation timestamp: {created_at!r}") from exc

    def _path(self, record):
        return self.directory / self._bucket(record.created_at) / f"{record.id}.json"

    def _run_lock_path(self, record):
        path = self._path(record)
        return path.with_name(f".{record.id}.lock")

    def _path_for_id(self, run_id):
        if not self.directory.is_dir():
            raise FileNotFoundError(run_id)
        filename = f"{run_id}.json"
        for day in sorted(self.directory.iterdir(), reverse=True):
            if day.is_dir():
                path = day / filename
                if path.is_file():
                    return path
        raise FileNotFoundError(run_id)

    @staticmethod
    def _decode(path):
        try:
            with path.open() as file:
                data = json.load(file)
            version = data.get("schema_version")
            if version != _RUN_RECORD_SCHEMA_VERSION:
                raise ValueError(
                    f"unsupported schema version {version!r}; "
                    f"expected {_RUN_RECORD_SCHEMA_VERSION}"
                )
            return RunRecord(**data)
        except CorruptRunRecord:
            raise
        except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
            raise CorruptRunRecord(path, str(exc)) from exc

    def _warn_corrupt(self, exc):
        key = os.fspath(exc.path)
        if key in self._warned_corrupt:
            return
        self._warned_corrupt.add(key)
        warnings.warn(str(exc), RuntimeWarning, stacklevel=3)

    @staticmethod
    def _json_default(value):
        if isinstance(value, os.PathLike):
            return os.fspath(value)
        raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")

    def _write_unlocked(self, record):
        durable_mkdir(self.directory, parents=True, exist_ok=True, mode=0o700)
        path = self._path(record)
        durable_mkdir(path.parent, parents=True, exist_ok=True, mode=0o700)
        durable_write_json(
            path,
            asdict(record),
            mode=0o600,
            sort_keys=True,
            default=self._json_default,
        )
        return record

    @contextmanager
    def submission_lock(self):
        durable_mkdir(self.directory, parents=True, exist_ok=True, mode=0o700)
        with self._lock, FileLock(self.directory / ".submit.lock"):
            yield

    def save(self, record: RunRecord):
        with self._lock:
            return self._write_unlocked(record)

    def update(self, record: RunRecord, updater):
        with self._lock, FileLock(self._run_lock_path(record)):
            path = self._path(record)
            latest = self._decode(path) if path.is_file() else record
            updated = updater(latest)
            if updated.id != latest.id or updated.created_at != latest.created_at:
                raise ValueError("Run update cannot change persistent identity")
            return self._write_unlocked(updated)

    def get(self, run_id):
        return self._decode(self._path_for_id(run_id))

    def _day_directories(self, *, newest):
        if not self.directory.is_dir():
            return ()
        days = [path for path in self.directory.iterdir() if path.is_dir()]
        return tuple(sorted(days, key=lambda path: path.name, reverse=newest))

    def _records_in_directory(self, directory, *, newest):
        records = []
        for path in directory.glob("*.json"):
            try:
                records.append(self._decode(path))
            except CorruptRunRecord as exc:
                self._warn_corrupt(exc)
        return sorted(records, key=lambda record: record.created_at, reverse=newest)

    def runs(self, *, active=None, limit=None, newest=False):
        if limit is not None and limit < 0:
            raise ValueError("limit must be non-negative")
        if limit == 0:
            return ()

        with self._lock:
            records = []
            for day in self._day_directories(newest=newest):
                for record in self._records_in_directory(day, newest=newest):
                    if active is True and record.state in _FINAL_STATES:
                        continue
                    if active is False and record.state not in _FINAL_STATES:
                        continue
                    records.append(record)
                    if limit is not None and len(records) >= limit:
                        return tuple(records)
            return tuple(records)


class ExecutionManager:
    """Persist public run lifecycle without task-specific knowledge."""

    def __init__(self, store=None, *, executor_pool=None, submitting_grace=60.0):
        if submitting_grace < 0:
            raise ValueError("submitting_grace must be non-negative")
        self.store = store if store is not None else RunStore()
        self.executor_pool = (
            executor_pool if executor_pool is not None else default_executor_pool()
        )
        self.submitting_grace = float(submitting_grace)

    def register_executor(self, executor: Executor, endpoint=None):
        endpoint = endpoint or endpoint_for_executor(executor)
        return self.executor_pool.get(endpoint, lambda: executor)

    def executor_for(self, record):
        if record.endpoint is None:
            raise KeyError(f"Run {record.id} has no persisted execution endpoint")
        return self.executor_pool.for_endpoint(record.endpoint)

    @staticmethod
    def _result_files_for_spec(spec):
        directory = os.path.abspath(os.fspath(spec.directory))
        input_file = _result_file_path(getattr(spec, "input_file", None), directory)
        output_file = _result_file_path(getattr(spec, "output_file", None), directory)
        input_files = _unique(
            (
                input_file,
                *(
                    _result_file_path(value, directory)
                    for value in getattr(spec, "input_files", ())
                ),
            )
        )
        output_files = _unique(
            (
                output_file,
                *(
                    _result_file_path(value, directory)
                    for value in getattr(spec, "output_files", ())
                ),
            )
        )
        return ResultFilesRecord(
            input_file=input_file,
            output_file=output_file,
            input_files=input_files,
            output_files=output_files,
            availability_state=ArtifactAvailabilityState.PENDING,
        )

    @staticmethod
    def _concrete_result_files(manifests):
        paths = []
        seen = set()
        for manifest in manifests:
            matches = (
                glob.glob(manifest, recursive=True)
                if glob.has_magic(manifest)
                else (manifest,)
            )
            for path in sorted(matches):
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
            if glob.has_magic(manifest):
                if not glob.glob(manifest, recursive=True):
                    missing.append(manifest)
            elif not os.path.isfile(manifest):
                missing.append(manifest)
        return _unique(missing)

    def _capture_result_files(self, result_files):
        if result_files is None or result_files.captured_at is not None:
            return result_files
        return replace(
            result_files,
            input_versions=self._result_files_snapshot(result_files.input_files),
            output_versions=self._result_files_snapshot(
                result_files.output_files, sync=True
            ),
            captured_at=_now(),
        )

    def input_file_is_current(self, run):
        record = self.store.get(run) if isinstance(run, str) else self.store.get(run.id)
        result_files = record.result_files
        return self._file_is_current(
            result_files,
            None if result_files is None else result_files.input_file,
            () if result_files is None else result_files.input_versions,
        )

    def output_file_is_current(self, run):
        record = self.store.get(run) if isinstance(run, str) else self.store.get(run.id)
        result_files = record.result_files
        return self._file_is_current(
            result_files,
            None if result_files is None else result_files.output_file,
            () if result_files is None else result_files.output_versions,
        )

    def new_record(
        self,
        *,
        task,
        directory,
        dataset=None,
        endpoint=None,
        submission=None,
        result_files=None,
        input_definition=None,
        run_id=None,
    ):
        now = _now()
        return RunRecord(
            id=run_id or uuid4().hex,
            state=RunState.NEW,
            created_at=now,
            updated_at=now,
            finished_at=None,
            task=str(task).upper(),
            dataset=None if dataset is None else str(dataset),
            directory=os.path.abspath(os.fspath(directory)),
            endpoint=endpoint,
            submission=submission,
            result_files=result_files,
            input_definition=input_definition,
        )

    @staticmethod
    def _merged_state(latest, requested):
        if requested is None:
            return latest.state
        requested = RunState(requested)
        if latest.state in _FINAL_STATES:
            return latest.state
        if requested == RunState.UNKNOWN or latest.state == RunState.UNKNOWN:
            return requested
        if _STATE_RANK[requested] < _STATE_RANK[latest.state]:
            return latest.state
        return requested

    @staticmethod
    def _merged_result_files(latest, requested):
        if requested is _UNCHANGED or requested is None:
            return latest.result_files
        if isinstance(requested, dict):
            requested = ResultFilesRecord(**requested)
        current = latest.result_files
        if current is None:
            return requested
        if current.captured_at is not None and requested.captured_at is not None:
            return current
        return requested

    def transition(
        self,
        record,
        *,
        state=None,
        reference=_UNCHANGED,
        return_code=_UNCHANGED,
        message=_UNCHANGED,
        failure_phase=_UNCHANGED,
        result_files=_UNCHANGED,
    ):
        def merge(latest):
            final_locked = latest.state in _FINAL_STATES
            now = _now()
            merged_state = self._merged_state(latest, state)
            changes = {"state": merged_state, "updated_at": now}
            if latest.state not in _FINAL_STATES and merged_state in _FINAL_STATES:
                changes["finished_at"] = now
            merged_result_files = self._merged_result_files(latest, result_files)
            if merged_result_files is not latest.result_files:
                changes["result_files"] = merged_result_files
            if (
                reference is not _UNCHANGED
                and reference is not None
                and latest.reference is None
            ):
                changes["reference"] = reference
            if return_code is not _UNCHANGED and return_code is not None:
                if not final_locked or latest.return_code is None:
                    changes["return_code"] = return_code
            if message is not _UNCHANGED and message is not None:
                if not final_locked or latest.message is None:
                    changes["message"] = message
            if failure_phase is not _UNCHANGED:
                normalized = _failure_phase(failure_phase)
                if not final_locked or latest.failure_phase is None:
                    changes["failure_phase"] = normalized
            return replace(latest, **changes)

        return self.store.update(record, merge)

    def update_artifact_availability(self, record, state, message=None):
        state = ArtifactAvailabilityState(state)

        def update(latest):
            result_files = latest.result_files
            if result_files is None:
                return latest
            return replace(
                latest,
                updated_at=_now(),
                result_files=replace(
                    result_files,
                    availability_state=state,
                    availability_message=message,
                ),
            )

        return self.store.update(record, update)

    def _is_stale_submitting(self, record, *, now=None):
        if record.state != RunState.SUBMITTING or record.reference is not None:
            return False
        now = now or datetime.now(timezone.utc)
        return (
            now - _timestamp(record.updated_at)
        ).total_seconds() >= self.submitting_grace

    def recover_submitting(self, record):
        record = self.store.get(record.id)
        if not self._is_stale_submitting(record):
            return record
        executor = self.executor_for(record)
        try:
            status = executor.recover_submission(
                record.id,
                record.submission,
                submitted_at=record.updated_at,
            )
        except NotImplementedError:
            return record
        if status is not None:
            return self.update(record, status)

        def mark_absent(latest):
            if latest.state != RunState.SUBMITTING or latest.reference is not None:
                return latest
            now = _now()
            return replace(
                latest,
                state=RunState.FAILED,
                updated_at=now,
                finished_at=now,
                message="Stale submission was not found by its execution backend",
                failure_phase=FailurePhase.SUBMISSION,
            )

        return self.store.update(record, mark_absent)

    def reserve(self, record):
        if record.submission is None:
            raise ValueError("Persistent run lacks a submission key")

        attempted_recovery = set()
        while True:
            with self.store.submission_lock():
                conflicts = tuple(
                    active
                    for active in self.store.runs(active=True)
                    if active.id != record.id
                    and active.submission is not None
                    and record.submission.conflicts_with(active.submission)
                )
                if not conflicts:
                    if record.result_files is not None:
                        self._sync_result_files(record.result_files.input_files)
                    return self.transition(
                        record,
                        state=RunState.SUBMITTING,
                        reference=None,
                    )
                stale = tuple(
                    conflict
                    for conflict in conflicts
                    if conflict.id not in attempted_recovery
                    and self._is_stale_submitting(conflict)
                )

            if not stale:
                raise DuplicateSubmissionError(record, conflicts)
            progressed = False
            for conflict in stale:
                attempted_recovery.add(conflict.id)
                recovered = self.recover_submitting(conflict)
                if recovered.state != RunState.SUBMITTING or recovered.reference is not None:
                    progressed = True
            if not progressed:
                raise DuplicateSubmissionError(record, conflicts)

    def update(self, record, status: RunStatus):
        phase = _failure_phase(status.failure_phase)
        return_code = status.return_code
        state = RunState.POSTPROCESSING if status.state == RunState.FINISHED else status.state
        if state == RunState.FAILED and return_code is None:
            return_code = 1
        return self.transition(
            record,
            state=state,
            reference=status.reference,
            return_code=return_code,
            message=status.message,
            failure_phase=phase,
        )

    def complete(
        self,
        record,
        *,
        state,
        reference=None,
        return_code=None,
        message=None,
        failure_phase=None,
        artifacts_available=False,
    ):
        state = RunState(state)
        phase = _failure_phase(failure_phase)
        result_files = record.result_files

        if artifacts_available and result_files is not None:
            result_files = replace(
                result_files,
                availability_state=ArtifactAvailabilityState.AVAILABLE,
                availability_message=None,
            )
            try:
                result_files = self._capture_result_files(result_files)
            except OSError as exc:
                state = RunState.FAILED
                phase = FailurePhase.POSTPROCESSING
                message = str(exc)
                result_files = replace(
                    result_files,
                    availability_state=ArtifactAvailabilityState.FAILED,
                    availability_message=str(exc),
                )
        elif phase in (FailurePhase.ARTIFACTS, FailurePhase.DOWNLOAD) and result_files is not None:
            result_files = replace(
                result_files,
                availability_state=ArtifactAvailabilityState.FAILED,
                availability_message=message,
            )

        return self.transition(
            record,
            state=state,
            reference=reference,
            return_code=return_code,
            message=message,
            failure_phase=phase,
            result_files=result_files,
        )

    def check_result_files(self, run):
        return self.store.get(run) if isinstance(run, str) else self.store.get(run.id)

    def retrieve_result_files(self, run):
        record = self.store.get(run) if isinstance(run, str) else self.store.get(run.id)
        if record.state not in _FINAL_STATES:
            raise RuntimeError(f"Run {record.id} is not finished")

        missing = self._missing_result_files(record.result_files)
        if not missing or record.reference is None:
            return record

        executor = self.executor_for(record)
        if not executor.supports_result_file_retrieval:
            return record

        executor.fetch_result_files(record.reference, missing)
        refreshed = self.store.get(record.id)
        remaining = self._missing_result_files(refreshed.result_files)
        if not remaining:
            return self.update_artifact_availability(
                refreshed, ArtifactAvailabilityState.AVAILABLE
            )
        return self.update_artifact_availability(
            refreshed,
            ArtifactAvailabilityState.FAILED,
            "Captured artifacts are still missing: " + ", ".join(remaining),
        )

    def _submission_context(self, executor, spec):
        provider = getattr(executor, "submission_context", None)
        return None if provider is None else provider(spec)

    @staticmethod
    def _complete(persistent, result, callback):
        persistent._sync_completion(result)
        if callback is not None:
            callback(result)

    def create_run(
        self,
        executor,
        spec,
        *,
        process_result,
        task,
        input_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
    ):
        """Create a persistent run from raw execution plus generic processing."""
        endpoint = endpoint_for_executor(executor)
        submission = SubmissionKey(
            endpoint=endpoint,
            directory=executor.submission_directory(spec),
            task=task,
            dataset=dataset,
            context=self._submission_context(executor, spec),
        )
        record = self.new_record(
            task=task,
            input_definition=input_definition,
            dataset=dataset,
            directory=spec.directory,
            endpoint=endpoint,
            submission=submission,
            result_files=self._result_files_for_spec(spec),
            run_id=run_id,
        )
        self.register_executor(executor, endpoint)

        raw = executor.create_persistent_run(spec, run_id=record.id)
        processed = PostprocessingRunner(raw, process_result)
        persistent = PersistentRunner(self, executor, processed, record)
        processed.callback = lambda result: self._complete(
            persistent, result, callback
        )
        return persistent

    def refresh(self, records=None):
        records = tuple(self.store.runs(active=True) if records is None else records)
        current = {record.id: record for record in records}
        groups = {}
        for record in records:
            if record.reference is None:
                continue
            executor = self.executor_for(record)
            groups.setdefault(executor, []).append(record)

        for executor, group in groups.items():
            statuses = executor.status(*(record.reference for record in group))
            if len(statuses) != len(group):
                raise RuntimeError("Executor returned wrong number of run statuses")
            for record, status in zip(group, statuses):
                current[record.id] = self.update(record, status)
        return tuple(current[record.id] for record in records)

    def attach(
        self,
        run,
        *,
        process_result,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        """Attach raw execution, then resume client-side postprocessing."""
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if executor is None:
            executor = self.executor_for(record)
        else:
            endpoint = record.endpoint or endpoint_for_executor(executor)
            self.register_executor(executor, endpoint)

        raw = executor.attach(
            record.reference,
            output_listener=output_listener,
            start=False,
        )
        processed = PostprocessingRunner(raw, process_result)
        persistent = PersistentRunner(self, executor, processed, record)
        processed.callback = lambda result: self._complete(
            persistent, result, callback
        )
        persistent.start()
        return persistent


class PersistentRunner:
    """Caller-owned runner that persists the public lifecycle directly."""

    def __init__(self, manager, executor, runner, record):
        self.manager = manager
        self.executor = executor
        self._runner = runner
        self._record = record

    @property
    def id(self):
        return self._record.id

    @property
    def record(self):
        return self._record

    @property
    def result(self):
        return self._runner.result

    @property
    def reference(self):
        return self._runner.reference

    @property
    def exception(self):
        return self._runner.exception

    @property
    def failure_phase(self):
        return _failure_phase(getattr(self._runner, "failure_phase", None))

    @property
    def state(self):
        state = self._runner.state
        if state == RunState.NEW and self._record.state != RunState.NEW:
            return self._record.state
        if self._record.state != RunState.NEW and state != self._record.state:
            self._record = self.manager.transition(
                self._record,
                state=state,
                reference=self.reference,
                failure_phase=self.failure_phase,
            )
        return self._record.state if self._record.state != RunState.NEW else state

    def add_output_listener(self, callback):
        return self._runner.add_output_listener(callback)

    def add_output_error_listener(self, callback):
        return self._runner.add_output_error_listener(callback)

    def submit(self):
        if self._record.state != RunState.NEW:
            raise RuntimeError(f"Run is already {self._record.state.value}")
        self._record = self.manager.reserve(self._record)
        try:
            result = self._runner.submit()
        except SubmissionNotStarted as exc:
            phase = self.failure_phase or FailurePhase.SUBMISSION
            message = getattr(self._runner, "failure_message", None) or str(exc)
            self._record = self.manager.complete(
                self._record,
                state=RunState.FAILED,
                message=message,
                failure_phase=phase,
                artifacts_available=False,
            )
            raise
        except BaseException:
            # No reference means submission may be uncertain. Keep SUBMITTING so
            # recover_submitting() can resolve it instead of inventing failure.
            reference = self.reference
            if reference is not None:
                self._record = self.manager.transition(
                    self._record,
                    state=self._runner.state,
                    reference=reference,
                    failure_phase=self.failure_phase,
                )
            raise

        reference = self.reference
        if reference is None:
            raise RuntimeError("Executor submit returned without a run reference")
        self.manager.register_executor(self.executor, self._record.endpoint)
        self._record = self.manager.transition(
            self._record,
            state=self._runner.state,
            reference=reference,
            failure_phase=self.failure_phase,
        )
        return result

    def start(self):
        self._runner.start()
        if self.reference is not None:
            self._record = self.manager.transition(
                self._record,
                state=self._runner.state,
                reference=self.reference,
                failure_phase=self.failure_phase,
            )

    def _sync_completion(self, result):
        state = self._runner.state
        phase = self.failure_phase
        message = getattr(self._runner, "failure_message", None)
        error = self._runner.exception
        if error is not None:
            message = str(error)
        if state not in _FINAL_STATES:
            return
        if state == RunState.FAILED and phase is None:
            phase = FailurePhase.EXECUTION

        artifacts_available = phase not in _NO_ARTIFACTS_FAILURES
        self._record = self.manager.complete(
            self._record,
            state=state,
            reference=self.reference,
            return_code=getattr(self._runner, "return_code", None),
            message=message,
            failure_phase=phase,
            artifacts_available=artifacts_available,
        )

    def wait(self):
        try:
            return self._runner.wait()
        finally:
            if self._record.state != RunState.NEW:
                self._sync_completion(self.result)

    def run(self):
        if self._record.state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        return await asyncio.to_thread(self.run)

    def refresh(self):
        self._record = self.manager.refresh((self._record,))[0]
        return self._record

    def cancel(self):
        try:
            return self._runner.cancel()
        finally:
            if self._record.state != RunState.NEW:
                self._record = self.manager.transition(
                    self._record,
                    state=self._runner.state,
                    reference=self.reference,
                    failure_phase=self.failure_phase,
                )
