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
from typing import Optional
from uuid import uuid4

from platformdirs import user_state_dir

from ..common.file_locks import FileLock
from ..common.file_versions import FileState, FileVersion
from .base import (
    Executor,
    RunReference,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .endpoint import ExecutorEndpoint, endpoint_for_executor
from .pool import default_executor_pool
from .result_capture import PersistentOutputCapture


_RUN_RECORD_SCHEMA_VERSION = 1
_FINAL_STATES = {
    RunState.AVAILABLE,
    RunState.CHECKED,
    RunState.FINISHED,
    RunState.EXECUTION_FAILED,
    RunState.FINALIZATION_FAILED,
    RunState.FAILED,
    RunState.CANCELLED,
}
_STATE_RANK = {
    RunState.NEW: 0,
    RunState.SUBMITTING: 1,
    RunState.QUEUED: 2,
    RunState.RUNNING: 3,
    RunState.FINALIZING: 4,
    RunState.AVAILABLE: 5,
    RunState.FINISHED: 5,
    RunState.CHECKED: 6,
    RunState.EXECUTION_FAILED: 5,
    RunState.FINALIZATION_FAILED: 5,
    RunState.FAILED: 5,
    RunState.CANCELLED: 5,
}
_UNCHANGED = object()


def _now():
    return datetime.now(timezone.utc).isoformat()


def _timestamp(value):
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed


def _result_file_path(value, directory):
    """Return a stable absolute local path for a declared result file."""
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


def _state_from_backend(state):
    """Map one backend status into the authoritative persistent lifecycle."""
    state = RunState(state)
    if state in (RunState.FINISHED, RunState.EXECUTION_FAILED, RunState.FAILED):
        return RunState.FINALIZING
    return state


def _completion_state(backend_state, error):
    """Resolve lifecycle state after execution/file finalization."""
    backend_state = RunState(backend_state)
    if backend_state in (
        RunState.FINISHED,
        RunState.EXECUTION_FAILED,
        RunState.FAILED,
    ):
        if error is not None:
            return RunState.FINALIZATION_FAILED
        if backend_state == RunState.FINISHED:
            return RunState.AVAILABLE
        return backend_state
    if backend_state == RunState.CANCELLED:
        return RunState.CANCELLED
    if error is not None:
        return RunState.FINALIZATION_FAILED
    return _state_from_backend(backend_state)


class CorruptRunRecord(RuntimeError):
    """A persisted run record cannot be decoded or validated."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Corrupt run record {self.path}: {message}")


class FileRetrievalState(str, Enum):
    """Availability of requested remote files in the local workspace."""

    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"


@dataclass(frozen=True)
class ResultFilesRecord:
    """Persistent manifests and exact file versions captured at finalization."""

    input_file: Optional[str] = None
    output_file: Optional[str] = None
    input_files: tuple[str, ...] = ()
    output_files: tuple[str, ...] = ()
    input_versions: tuple[FileVersion, ...] = ()
    output_versions: tuple[FileVersion, ...] = ()
    captured_at: Optional[str] = None
    retrieval_state: Optional[FileRetrievalState] = None
    retrieval_message: Optional[str] = None

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
        if self.retrieval_state is not None and not isinstance(
            self.retrieval_state, FileRetrievalState
        ):
            object.__setattr__(
                self,
                "retrieval_state",
                FileRetrievalState(self.retrieval_state),
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
    """Persistent namespace used only for duplicate-submit protection."""

    endpoint: ExecutorEndpoint
    directory: str
    task: str
    dataset: Optional[str] = None

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
    """Persistent identity and last known lifecycle state of one run.

    Terminal backend completion first becomes ``FINALIZING``. Successful file
    finalization becomes ``AVAILABLE``; parsing/checking then advances to
    ``CHECKED``. Retryable check failures stay ``AVAILABLE`` with a failure
    reason, while unexpected check failures become ``FAILED``. Failed/cancelled
    executions are never checked automatically.
    """

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
        if self.schema_version != _RUN_RECORD_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported schema version {self.schema_version}; "
                f"expected {_RUN_RECORD_SCHEMA_VERSION}"
            )

    @property
    def failure_reason(self):
        """Current lifecycle/check failure reason without duplicate persistence."""
        return self.message


class RunStore:
    """Persistent run history stored as date-sharded atomic JSON files."""

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
            data.pop("requested_task", None)
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
        self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        path = self._path(record)
        path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
        try:
            temporary.write_text(
                json.dumps(asdict(record), sort_keys=True, default=self._json_default)
            )
            try:
                os.chmod(temporary, 0o600)
            except OSError:
                pass
            os.replace(temporary, path)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        return record

    @contextmanager
    def submission_lock(self):
        self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
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
    """Orchestrate optional persistent run bookkeeping."""

    def __init__(self, store=None, *, executor_pool=None, submitting_grace=60.0):
        if submitting_grace < 0:
            raise ValueError("submitting_grace must be non-negative")
        self.store = store if store is not None else RunStore()
        self.executor_pool = executor_pool if executor_pool is not None else default_executor_pool()
        self.submitting_grace = float(submitting_grace)
        self._legacy_executors = {}
        self._lock = RLock()

    @staticmethod
    def _reference_key(reference):
        return reference.executor, reference.computer

    def register_executor(self, executor: Executor, endpoint=None, reference=None):
        endpoint = endpoint or endpoint_for_executor(executor)
        pooled = self.executor_pool.get(endpoint, lambda: executor)
        if reference is not None:
            with self._lock:
                self._legacy_executors[self._reference_key(reference)] = pooled
        return pooled

    def executor_for(self, record):
        if record.endpoint is not None:
            return self.executor_pool.for_endpoint(record.endpoint)
        reference = record.reference
        if reference is None:
            raise KeyError(f"Run {record.id} has neither endpoint nor backend reference")
        with self._lock:
            executor = self._legacy_executors.get(self._reference_key(reference))
            if executor is None and reference.computer is not None:
                executor = self._legacy_executors.get((reference.executor, None))
        if executor is None:
            raise KeyError(
                f"Run {record.id} predates persisted endpoint identity and no live executor is registered for it"
            )
        return executor

    @staticmethod
    def _result_files_for_spec(spec):
        directory = os.path.abspath(os.fspath(spec.directory))
        input_file = _result_file_path(getattr(spec, "input_file", None), directory)
        output_file = _result_file_path(getattr(spec, "output_file", None), directory)
        input_files = _unique((input_file, *(_result_file_path(value, directory) for value in getattr(spec, "input_files", ()))))
        output_files = _unique((output_file, *(_result_file_path(value, directory) for value in getattr(spec, "output_files", ()))))
        return ResultFilesRecord(
            input_file=input_file,
            output_file=output_file,
            input_files=input_files,
            output_files=output_files,
        )

    @staticmethod
    def _result_files_snapshot(manifests):
        versions = []
        seen = set()
        for manifest in manifests:
            matches = glob.glob(manifest, recursive=True) if glob.has_magic(manifest) else (manifest,)
            for path in sorted(matches):
                path = os.path.abspath(path)
                if path in seen or not os.path.isfile(path):
                    continue
                seen.add(path)
                versions.append(FileVersion.capture(path, require_digest=True))
        return tuple(sorted(versions, key=lambda value: os.fspath(value.path)))

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

    @classmethod
    def _result_files_missing(cls, result_files):
        return bool(cls._missing_result_files(result_files))

    def _capture_result_files(self, result_files):
        if result_files is None or result_files.captured_at is not None:
            return result_files
        return replace(
            result_files,
            input_versions=self._result_files_snapshot(result_files.input_files),
            output_versions=self._result_files_snapshot(result_files.output_files),
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
            if reference is not _UNCHANGED and reference is not None and latest.reference is None:
                changes["reference"] = reference
            if return_code is not _UNCHANGED and return_code is not None:
                if not final_locked or latest.return_code is None:
                    changes["return_code"] = return_code
            if message is not _UNCHANGED and message is not None:
                if not final_locked or latest.message is None:
                    changes["message"] = message
            return replace(latest, **changes)
        return self.store.update(record, merge)

    def fail_finalization(self, record, error):
        record = self.store.get(record.id)
        if record.state != RunState.FINALIZING:
            return record
        return self.finalize(record, state=RunState.FINALIZATION_FAILED, message=str(error))

    def retry_finalization(self, record):
        def retry(latest):
            if latest.state != RunState.FINALIZATION_FAILED:
                raise RuntimeError(
                    f"Run {latest.id} is {latest.state.value}, not a finalization failure"
                )
            return replace(
                latest,
                state=RunState.FINALIZING,
                updated_at=_now(),
                finished_at=None,
                message=None,
            )
        return self.store.update(record, retry)

    def mark_available(self, record):
        def update(latest):
            if latest.state == RunState.AVAILABLE:
                return latest
            if latest.state != RunState.FINISHED:
                raise RuntimeError(
                    f"Run {latest.id} is {latest.state.value}, not legacy finished"
                )
            return replace(latest, state=RunState.AVAILABLE, updated_at=_now())
        return self.store.update(record, update)

    def mark_checked(self, record):
        def update(latest):
            if latest.state == RunState.CHECKED and latest.message is None:
                return latest
            if latest.state not in (RunState.AVAILABLE, RunState.CHECKED):
                raise RuntimeError(f"Run {latest.id} is {latest.state.value}, not checkable")
            return replace(
                latest,
                state=RunState.CHECKED,
                updated_at=_now(),
                message=None,
            )
        return self.store.update(record, update)

    def record_result_check_failure(self, record, error, *, retryable):
        def update(latest):
            if latest.state not in (RunState.AVAILABLE, RunState.CHECKED):
                return latest
            return replace(
                latest,
                state=RunState.AVAILABLE if retryable else RunState.FAILED,
                updated_at=_now(),
                message=str(error),
            )
        return self.store.update(record, update)

    def _is_stale_submitting(self, record, *, now=None):
        if record.state != RunState.SUBMITTING or record.reference is not None:
            return False
        now = now or datetime.now(timezone.utc)
        return (now - _timestamp(record.updated_at)).total_seconds() >= self.submitting_grace

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
                    return self.transition(record, state=RunState.SUBMITTING, reference=None)
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
        backend_state = RunState(status.state)
        return_code = status.return_code
        if backend_state in (RunState.EXECUTION_FAILED, RunState.FAILED) and return_code is None:
            return_code = 1
        return self.transition(
            record,
            state=_state_from_backend(backend_state),
            reference=status.reference,
            return_code=return_code,
            message=status.message,
        )

    def finalize(
        self,
        record,
        *,
        state,
        reference=None,
        return_code=None,
        message=None,
    ):
        state = RunState(state)
        result_files = record.result_files
        if state != RunState.FINALIZATION_FAILED:
            result_files = self._capture_result_files(result_files)
        return self.transition(
            record,
            state=state,
            reference=reference,
            return_code=return_code,
            message=message,
            result_files=result_files,
        )

    def update_file_retrieval(self, record, state, message=None):
        state = FileRetrievalState(state)
        def update(latest):
            result_files = latest.result_files
            if result_files is None:
                return latest
            return replace(
                latest,
                updated_at=_now(),
                result_files=replace(
                    result_files,
                    retrieval_state=state,
                    retrieval_message=message,
                ),
            )
        return self.store.update(record, update)

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
            return self.update_file_retrieval(refreshed, FileRetrievalState.AVAILABLE)
        if refreshed.state in (RunState.AVAILABLE, RunState.CHECKED, RunState.FINISHED):
            return self.update_file_retrieval(
                refreshed,
                FileRetrievalState.FAILED,
                "Captured remote files are still missing: " + ", ".join(remaining),
            )
        return refreshed

    @staticmethod
    def _complete(persistent, result, callback):
        first_available = persistent._sync_completion(result)
        delivered = result
        if first_available:
            delivered = persistent._check_result()
        if callback is not None:
            callback(delivered)

    def _wrap_runner(self, executor, record, reader, factory, callback):
        guard = RLock()
        holder = {"runner": None, "pending": _UNCHANGED}

        def completed(result):
            with guard:
                persistent = holder["runner"]
                if persistent is None:
                    holder["pending"] = result
                    return
            self._complete(persistent, result, callback)

        runner = factory(completed)
        persistent = PersistentRunner(self, executor, runner, record, reader)
        with guard:
            holder["runner"] = persistent
            pending = holder["pending"]
            holder["pending"] = _UNCHANGED
        if pending is not _UNCHANGED:
            self._complete(persistent, pending, callback)
        return persistent

    def create_run(
        self,
        executor,
        reader,
        spec,
        *,
        task,
        input_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
    ):
        endpoint = endpoint_for_executor(executor)
        submission = SubmissionKey(
            endpoint=endpoint,
            directory=executor.submission_directory(spec),
            task=task,
            dataset=dataset,
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
        capture = PersistentOutputCapture(reader)
        return self._wrap_runner(
            executor,
            record,
            reader,
            lambda completed: executor.create_persistent_run(
                capture,
                spec,
                run_id=record.id,
                callback=completed,
            ),
            callback,
        )

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
        reader,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if executor is None:
            executor = self.executor_for(record)
        else:
            endpoint = record.endpoint or endpoint_for_executor(executor)
            self.register_executor(executor, endpoint, record.reference)
        capture = PersistentOutputCapture(reader)
        return self._wrap_runner(
            executor,
            record,
            reader,
            lambda completed: executor.attach(
                record.reference,
                capture,
                callback=completed,
                output_listener=output_listener,
            ),
            callback,
        )


class PersistentRunner:
    """Caller-owned runner that persists execution and result lifecycle transitions."""

    def __init__(self, manager, executor, runner, record, reader):
        self.manager = manager
        self.executor = executor
        self._runner = runner
        self._record = record
        self._reader = reader
        self._result = runner.result

    @property
    def id(self):
        return self._record.id

    @property
    def record(self):
        return self._record

    @property
    def result(self):
        return self._result

    @property
    def reference(self):
        return self._runner.reference

    @property
    def state(self):
        if self._record.state in _FINAL_STATES:
            return self._record.state
        backend_state = self._runner.state
        state = _state_from_backend(backend_state)
        if backend_state == RunState.NEW and self._record.state != RunState.NEW:
            state = self._record.state
        if self._record.state != RunState.NEW and state != self._record.state:
            self._record = self.manager.transition(
                self._record,
                state=state,
                reference=self.reference,
            )
            state = self._record.state
        return state

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
            self._record = self.manager.transition(
                self._record,
                state=RunState.FAILED,
                message=str(exc),
            )
            raise
        except BaseException:
            reference = self.reference
            if reference is not None:
                backend_state = self._runner.state
                if backend_state == RunState.NEW:
                    backend_state = RunState.QUEUED
                self._record = self.manager.transition(
                    self._record,
                    state=_state_from_backend(backend_state),
                    reference=reference,
                )
            raise
        reference = self.reference
        if reference is None:
            raise RuntimeError("Executor submit returned without a run reference")
        self.manager.register_executor(self.executor, self._record.endpoint, reference)
        backend_state = self._runner.state
        if backend_state == RunState.NEW:
            backend_state = RunState.QUEUED
        self._record = self.manager.transition(
            self._record,
            state=_state_from_backend(backend_state),
            reference=reference,
        )
        self._sync_file_retrieval()
        return result

    def _sync_file_retrieval(self):
        state = getattr(self._runner, "file_retrieval_state", None)
        if state is None:
            return
        self._record = self.manager.update_file_retrieval(
            self._record,
            state,
            getattr(self._runner, "file_retrieval_message", None),
        )

    def _sync_completion(self, result):
        previous = self.manager.store.get(self._record.id)
        self._sync_file_retrieval()
        backend_state = self._runner.state
        error = self._runner.exception
        state = _completion_state(backend_state, error)
        return_code = getattr(result, "return_code", None)
        if return_code is None:
            return_code = getattr(self._runner, "return_code", None)
        self._record = self.manager.finalize(
            self._record,
            state=state,
            reference=self.reference,
            return_code=return_code,
            message=None if error is None else str(error),
        )
        return previous.state not in _FINAL_STATES and self._record.state == RunState.AVAILABLE

    def _check_result(self):
        from ..outputs.persistent import check_persistent_result
        self._result = check_persistent_result(
            self.manager,
            self._record,
            reader=self._reader,
        )
        self._record = self.manager.store.get(self._record.id)
        return self._result

    def wait(self):
        try:
            result = self._runner.wait()
            return self._result if self._record.state == RunState.CHECKED else result
        finally:
            if self._record.state != RunState.NEW and self._record.state not in _FINAL_STATES:
                self._sync_completion(self._runner.result)

    def run(self):
        if self._record.state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.run)

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
                    state=_state_from_backend(self._runner.state),
                    reference=self.reference,
                )
