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
from .base import (
    Executor,
    RunReference,
    RunState,
    RunStatus,
    SubmissionNotStarted,
)
from .endpoint import ExecutorEndpoint, endpoint_for_executor
from .pool import default_executor_pool


_RUN_RECORD_SCHEMA_VERSION = 2
_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}
_STATE_RANK = {
    RunState.NEW: 0,
    RunState.SUBMITTING: 1,
    RunState.QUEUED: 2,
    RunState.RUNNING: 3,
    RunState.FINISHING: 4,
    RunState.FINISHED: 5,
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


def _processing_state(backend_state):
    """Map a terminal backend outcome to the still-running finalization phase."""
    backend_state = RunState(backend_state)
    if backend_state in (RunState.FINISHED, RunState.FAILED):
        return RunState.FINISHING
    return backend_state


class CorruptRunRecord(RuntimeError):
    """A persisted run record cannot be decoded or validated."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Corrupt run record {self.path}: {message}")


class ResultFilesState(str, Enum):
    """Freshness of files produced by a completed calculation."""

    CURRENT = "current"
    OVERWRITTEN = "overwritten"


class FileRetrievalState(str, Enum):
    """Availability of requested remote files in the local workspace."""

    PENDING = "pending"
    AVAILABLE = "available"
    FAILED = "failed"


@dataclass(frozen=True)
class FileFingerprint:
    """Cheap local identity used to notice later file changes."""

    path: str
    size: int
    mtime_ns: int

    def __post_init__(self):
        object.__setattr__(self, "path", os.fspath(self.path))


@dataclass(frozen=True)
class ResultFilesRecord:
    """Persistent file manifests and snapshots captured at finalization.

    ``fingerprints`` describe calculation outputs and determine ``state``.
    ``input_fingerprints`` preserve reconstructable historical input context,
    but an unavailable/changed input does not invalidate output that is still
    current. Cache files are intentionally absent because reusable caches may
    legitimately be updated by later calculations.
    """

    input_file: Optional[str] = None
    output_file: Optional[str] = None
    input_files: tuple[str, ...] = ()
    output_files: tuple[str, ...] = ()
    input_fingerprints: tuple[FileFingerprint, ...] = ()
    fingerprints: tuple[FileFingerprint, ...] = ()
    captured_at: Optional[str] = None
    state: ResultFilesState = ResultFilesState.CURRENT
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
        for name in ("input_fingerprints", "fingerprints"):
            object.__setattr__(
                self,
                name,
                tuple(
                    value
                    if isinstance(value, FileFingerprint)
                    else FileFingerprint(**value)
                    for value in getattr(self, name)
                ),
            )
        if not isinstance(self.state, ResultFilesState):
            object.__setattr__(self, "state", ResultFilesState(self.state))
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
    def all_fingerprints(self):
        return tuple(
            sorted(
                {
                    item.path: item
                    for item in (*self.input_fingerprints, *self.fingerprints)
                }.values(),
                key=lambda item: item.path,
            )
        )


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

    ``task`` is the concrete native execution task. ``task_definition`` is the
    bound ASE2SPRKKR task identity when it was known at submission time. Keeping
    both allows result reconstruction without conflating task policy with a
    concrete ``TASK.TASK`` submode.

    ``state`` is the public run lifecycle, including result finalization. A
    terminal backend outcome is therefore represented as ``FINISHING`` until
    output parsing and the completion callback have finished. A finalization
    error makes the public state ``FAILED`` even if the backend itself exited
    successfully. ``backend_state`` records that independent backend outcome.

    ``finished_at`` records the first terminal lifecycle transition and is
    immutable afterwards. Later metadata changes, including result-file
    freshness and remote file retrieval, update ``updated_at`` but never
    ``finished_at``.
    """

    id: str
    state: RunState
    created_at: str
    updated_at: str
    task: str
    directory: str
    task_definition: Optional[str] = None
    finished_at: Optional[str] = None
    dataset: Optional[str] = None
    endpoint: Optional[ExecutorEndpoint] = None
    submission: Optional[SubmissionKey] = None
    reference: Optional[RunReference] = None
    backend_state: Optional[RunState] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    result_files: Optional[ResultFilesRecord] = None
    schema_version: int = _RUN_RECORD_SCHEMA_VERSION

    def __post_init__(self):
        object.__setattr__(self, "task", str(self.task).upper())
        if self.task_definition is not None:
            object.__setattr__(
                self, "task_definition", str(self.task_definition).upper()
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
        if self.backend_state is not None and not isinstance(self.backend_state, RunState):
            object.__setattr__(self, "backend_state", RunState(self.backend_state))
        if self.schema_version != _RUN_RECORD_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported schema version {self.schema_version}; "
                f"expected {_RUN_RECORD_SCHEMA_VERSION}"
            )


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
    def _migrate_result_files(data):
        result_files = data.get("result_files")
        if result_files is None:
            return
        result_files.pop("output_digest", None)
        if "input_fingerprints" in result_files:
            return

        # Schema 1 stored inputs and outputs in one snapshot. Split exact input
        # manifests where possible; leave ambiguous paths in the output snapshot
        # so migration remains conservative rather than accepting changed output.
        fingerprints = result_files.get("fingerprints", ())
        input_exact = {
            os.path.abspath(path)
            for path in result_files.get("input_files", ())
            if not glob.has_magic(path)
        }
        output_exact = {
            os.path.abspath(path)
            for path in result_files.get("output_files", ())
            if not glob.has_magic(path)
        }
        input_fingerprints = []
        output_fingerprints = []
        for fingerprint in fingerprints:
            path = os.path.abspath(fingerprint["path"])
            if path in input_exact and path not in output_exact:
                input_fingerprints.append(fingerprint)
            else:
                output_fingerprints.append(fingerprint)
        result_files["input_fingerprints"] = input_fingerprints
        result_files["fingerprints"] = output_fingerprints

    @staticmethod
    def _decode(path):
        try:
            with path.open() as file:
                data = json.load(file)
            data.setdefault("finished_at", None)
            version = data.get("schema_version", 1)
            if version == 1:
                data["schema_version"] = _RUN_RECORD_SCHEMA_VERSION
                data.setdefault("task_definition", None)
            elif version == _RUN_RECORD_SCHEMA_VERSION:
                data.setdefault("task_definition", None)
            RunStore._migrate_result_files(data)
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

    def _write_unlocked(self, record):
        self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        path = self._path(record)
        path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
        try:
            temporary.write_text(json.dumps(asdict(record), sort_keys=True))
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
        """Serialize duplicate checking and the SUBMITTING record creation."""
        self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        with self._lock, FileLock(self.directory / ".submit.lock"):
            yield

    def save(self, record: RunRecord):
        """Atomically store a complete record supplied by the caller."""
        with self._lock:
            return self._write_unlocked(record)

    def update(self, record: RunRecord, updater):
        """Read latest, merge under a per-run lock, and atomically replace it."""
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
        """Return stored runs without querying backends or touching result files."""
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
    """Orchestrate optional persistent run bookkeeping.

    The manager owns persistent records, not live runners. Returned
    :class:`PersistentRunner` objects remain owned by their callers. Executors
    are process-wide pooled by persisted endpoint identity, independently of
    this manager and its store.
    """

    def __init__(self, store=None, *, executor_pool=None, submitting_grace=60.0):
        if submitting_grace < 0:
            raise ValueError("submitting_grace must be non-negative")
        self.store = store if store is not None else RunStore()
        self.executor_pool = (
            executor_pool if executor_pool is not None else default_executor_pool()
        )
        self.submitting_grace = float(submitting_grace)
        self._legacy_executors = {}
        self._lock = RLock()

    @staticmethod
    def _reference_key(reference):
        return reference.executor, reference.computer

    def register_executor(self, executor: Executor, endpoint=None, reference=None):
        """Register a live executor in the independent endpoint pool."""
        endpoint = endpoint or endpoint_for_executor(executor)
        pooled = self.executor_pool.get(endpoint, lambda: executor)
        if reference is not None:
            with self._lock:
                self._legacy_executors[self._reference_key(reference)] = pooled
        return pooled

    def executor_for(self, record):
        """Return or reconstruct the executor needed for a persisted run."""
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
                f"Run {record.id} predates persisted endpoint identity and no live "
                "executor is registered for it"
            )
        return executor

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
        )

    @staticmethod
    def _result_files_snapshot(manifests):
        fingerprints = []
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
                stat = os.stat(path)
                fingerprints.append(
                    FileFingerprint(
                        path=path,
                        size=stat.st_size,
                        mtime_ns=stat.st_mtime_ns,
                    )
                )
        return tuple(sorted(fingerprints, key=lambda value: value.path))

    @staticmethod
    def _fingerprint_current(fingerprint):
        try:
            stat = os.stat(fingerprint.path)
        except FileNotFoundError:
            return False
        return (
            stat.st_size == fingerprint.size
            and stat.st_mtime_ns == fingerprint.mtime_ns
        )

    @staticmethod
    def _missing_result_files(result_files):
        """Return the concrete historical files that should be restored."""
        if result_files is None:
            return ()
        missing = [
            item.path
            for item in result_files.all_fingerprints
            if not os.path.isfile(item.path)
        ]
        for manifest in result_files.manifests:
            if glob.has_magic(manifest):
                # Once a snapshot exists, concrete fingerprint paths are the
                # historical authority. Re-expanding a wildcard later could
                # pull files produced by another calculation.
                if result_files.captured_at is None and not glob.glob(
                    manifest, recursive=True
                ):
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
            input_fingerprints=self._result_files_snapshot(result_files.input_files),
            fingerprints=self._result_files_snapshot(result_files.output_files),
            captured_at=_now(),
            state=ResultFilesState.CURRENT,
        )

    def input_file_is_current(self, run):
        """Whether the primary historical input is safe to use for reconstruction."""
        record = (
            self.store.get(run)
            if isinstance(run, str)
            else self.store.get(run.id)
        )
        result_files = record.result_files
        if result_files is None or result_files.input_file is None:
            return False
        path = os.path.abspath(result_files.input_file)
        if result_files.captured_at is None:
            return os.path.isfile(path)
        for fingerprint in result_files.input_fingerprints:
            if os.path.abspath(fingerprint.path) == path:
                return self._fingerprint_current(fingerprint)
        return False

    def new_record(
        self,
        *,
        task,
        directory,
        task_definition=None,
        dataset=None,
        endpoint=None,
        submission=None,
        result_files=None,
        run_id=None,
    ):
        """Create the in-memory identity of an unsubmitted persistent run."""
        now = _now()
        return RunRecord(
            id=run_id or uuid4().hex,
            state=RunState.NEW,
            created_at=now,
            updated_at=now,
            finished_at=None,
            task=str(task).upper(),
            task_definition=(
                None if task_definition is None else str(task_definition).upper()
            ),
            dataset=None if dataset is None else str(dataset),
            directory=os.path.abspath(os.fspath(directory)),
            endpoint=endpoint,
            submission=submission,
            result_files=result_files,
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
    def _merged_backend_state(latest, requested):
        if requested is _UNCHANGED or requested is None:
            return latest.backend_state
        requested = RunState(requested)
        if latest.backend_state in _FINAL_STATES:
            return latest.backend_state
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
        if current.state == ResultFilesState.OVERWRITTEN:
            return current
        if current.captured_at is not None and requested.captured_at is not None:
            return current
        return requested

    def transition(
        self,
        record,
        *,
        state=None,
        reference=_UNCHANGED,
        backend_state=_UNCHANGED,
        return_code=_UNCHANGED,
        message=_UNCHANGED,
        result_files=_UNCHANGED,
    ):
        """Merge one lifecycle update into the latest persisted record."""

        def merge(latest):
            final_locked = latest.state in _FINAL_STATES
            now = _now()
            merged_state = self._merged_state(latest, state)
            changes = {
                "state": merged_state,
                "updated_at": now,
            }
            if latest.state not in _FINAL_STATES and merged_state in _FINAL_STATES:
                changes["finished_at"] = now
            merged_backend_state = self._merged_backend_state(latest, backend_state)
            if merged_backend_state is not latest.backend_state:
                changes["backend_state"] = merged_backend_state
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
            return replace(latest, **changes)

        return self.store.update(record, merge)

    def _is_stale_submitting(self, record, *, now=None):
        if record.state != RunState.SUBMITTING or record.reference is not None:
            return False
        now = now or datetime.now(timezone.utc)
        return (
            now - _timestamp(record.updated_at)
        ).total_seconds() >= self.submitting_grace

    def recover_submitting(self, record):
        """Recover one stale SUBMITTING record without guessing on uncertainty."""
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
        """Atomically reserve *record*'s persistent submission namespace."""
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
                if (
                    recovered.state != RunState.SUBMITTING
                    or recovered.reference is not None
                ):
                    progressed = True
            if not progressed:
                raise DuplicateSubmissionError(record, conflicts)

    def update(self, record, status: RunStatus):
        """Persist a backend status while keeping result finalization explicit."""
        backend_state = RunState(status.state)
        return self.transition(
            record,
            state=_processing_state(backend_state),
            backend_state=backend_state,
            reference=status.reference,
            return_code=status.return_code,
            message=status.message,
        )

    def finalize(
        self,
        record,
        *,
        state,
        reference=None,
        return_code=None,
        backend_state=None,
        message=None,
    ):
        """Capture result files while preserving the execution lifecycle outcome."""
        result_files = self._capture_result_files(record.result_files)
        if backend_state is None:
            backend_state = state if state in _FINAL_STATES else _UNCHANGED
        return self.transition(
            record,
            state=state,
            backend_state=backend_state,
            reference=reference,
            return_code=return_code,
            message=message,
            result_files=result_files,
        )

    def update_file_retrieval(self, record, state, message=None):
        """Persist local file availability without changing calculation state."""
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
        """Refresh output freshness without changing calculation outcome."""
        record = (
            self.store.get(run)
            if isinstance(run, str)
            else self.store.get(run.id)
        )

        def verify(latest):
            result_files = latest.result_files
            if (
                latest.state not in _FINAL_STATES
                or result_files is None
                or result_files.captured_at is None
                or result_files.state == ResultFilesState.OVERWRITTEN
            ):
                return latest
            if (
                self._result_files_snapshot(result_files.output_files)
                == result_files.fingerprints
            ):
                return latest
            return replace(
                latest,
                updated_at=_now(),
                result_files=replace(
                    result_files, state=ResultFilesState.OVERWRITTEN
                ),
            )

        return self.store.update(record, verify)

    def retrieve_result_files(self, run):
        """Restore locally missing files belonging to one terminal persistent run."""
        record = (
            self.store.get(run)
            if isinstance(run, str)
            else self.store.get(run.id)
        )
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
            return self.update_file_retrieval(
                refreshed, FileRetrievalState.AVAILABLE
            )
        if refreshed.state == RunState.FINISHED:
            return self.update_file_retrieval(
                refreshed,
                FileRetrievalState.FAILED,
                "Declared remote files are still missing: " + ", ".join(remaining),
            )
        return refreshed

    @staticmethod
    def _complete(persistent, result, callback):
        persistent._sync_completion(result)
        if callback is not None:
            callback(result)

    def _wrap_runner(self, executor, record, factory, callback):
        """Create a persistent wrapper without assuming callback timing."""
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
        persistent = PersistentRunner(self, executor, runner, record)
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
        task_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
    ):
        """Create a caller-owned runner whose lifecycle is persisted."""
        endpoint = endpoint_for_executor(executor)
        submission = SubmissionKey(
            endpoint=endpoint,
            directory=executor.submission_directory(spec),
            task=task,
            dataset=dataset,
        )
        record = self.new_record(
            task=task,
            task_definition=task_definition,
            dataset=dataset,
            directory=spec.directory,
            endpoint=endpoint,
            submission=submission,
            result_files=self._result_files_for_spec(spec),
            run_id=run_id,
        )
        self.register_executor(executor, endpoint)
        return self._wrap_runner(
            executor,
            record,
            lambda completed: executor.create_persistent_run(
                reader,
                spec,
                run_id=record.id,
                callback=completed,
            ),
            callback,
        )

    def refresh(self, records=None):
        """Refresh referenced runs, reconstructing endpoints when necessary."""
        records = tuple(
            self.store.runs(active=True) if records is None else records
        )
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
        """Attach to a stored run using an explicitly supplied output reader.

        ``executor`` is retained as an explicit override for records whose
        endpoint cannot be reconstructed from :class:`ExecutorEndpoint` alone.

        ``output_listener`` is installed by the backend before monitoring starts.
        Attached runners therefore cannot lose output replayed immediately during
        attach, including output from an already completed run.
        """
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if executor is None:
            executor = self.executor_for(record)
        else:
            endpoint = record.endpoint or endpoint_for_executor(executor)
            self.register_executor(executor, endpoint, record.reference)
        return self._wrap_runner(
            executor,
            record,
            lambda completed: executor.attach(
                record.reference,
                reader,
                callback=completed,
                output_listener=output_listener,
            ),
            callback,
        )


class PersistentRunner:
    """Caller-owned runner that persists execution lifecycle transitions."""

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
    def state(self):
        backend_state = self._runner.state
        state = _processing_state(backend_state)
        if backend_state == RunState.NEW and self._record.state != RunState.NEW:
            state = self._record.state
        if self._record.state != RunState.NEW and state != self._record.state:
            self._record = self.manager.transition(
                self._record,
                state=state,
                backend_state=backend_state,
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
                    state=_processing_state(backend_state),
                    backend_state=backend_state,
                    reference=reference,
                )
            raise

        reference = self.reference
        if reference is None:
            raise RuntimeError("Executor submit returned without a run reference")
        self.manager.register_executor(
            self.executor,
            self._record.endpoint,
            reference,
        )
        backend_state = self._runner.state
        if backend_state == RunState.NEW:
            backend_state = RunState.QUEUED
        self._record = self.manager.transition(
            self._record,
            state=_processing_state(backend_state),
            backend_state=backend_state,
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
        self._sync_file_retrieval()
        backend_state = self._runner.state
        error = self._runner.exception
        state = RunState.FAILED if error is not None else backend_state
        self._record = self.manager.finalize(
            self._record,
            state=state,
            backend_state=backend_state,
            reference=self.reference,
            return_code=getattr(result, "return_code", None),
            message=None if error is None else str(error),
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
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.run)

    def refresh(self):
        """Refresh this run's persisted backend status."""
        self._record = self.manager.refresh((self._record,))[0]
        return self._record

    def cancel(self):
        try:
            return self._runner.cancel()
        finally:
            if self._record.state != RunState.NEW:
                backend_state = self._runner.state
                self._record = self.manager.transition(
                    self._record,
                    state=_processing_state(backend_state),
                    backend_state=backend_state,
                    reference=self.reference,
                )
