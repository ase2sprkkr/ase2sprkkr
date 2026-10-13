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

from .base import Executor, RunReference, RunState, RunStatus
from .endpoint import ExecutorEndpoint, endpoint_for_executor
from .pool import default_executor_pool


_RUN_RECORD_SCHEMA_VERSION = 1
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


def _artifact_path(value, directory):
    """Return a stable absolute local path for a declared artifact."""
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


@contextmanager
def _filesystem_lock(path):
    """Hold an exclusive one-byte/file lock on *path* across processes."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    with path.open("a+b") as lock_file:
        if os.name == "nt":
            import msvcrt

            lock_file.seek(0, os.SEEK_END)
            if lock_file.tell() == 0:
                lock_file.write(b"\0")
                lock_file.flush()
            lock_file.seek(0)
            msvcrt.locking(lock_file.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                lock_file.seek(0)
                msvcrt.locking(lock_file.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


class CorruptRunRecord(RuntimeError):
    """A persisted run record cannot be decoded or validated."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Corrupt run record {self.path}: {message}")


class ArtifactState(str, Enum):
    """Freshness of files belonging to a completed calculation."""

    CURRENT = "current"
    OVERWRITTEN = "overwritten"


@dataclass(frozen=True)
class FileFingerprint:
    """Cheap local identity used to notice later artifact changes."""

    path: str
    size: int
    mtime_ns: int

    def __post_init__(self):
        object.__setattr__(self, "path", os.fspath(self.path))


@dataclass(frozen=True)
class ArtifactRecord:
    """Persistent artifact manifest and the snapshot captured at finalization.

    Cache files are intentionally absent. Reusable caches may legitimately be
    updated by later calculations and therefore never make a historical result
    ``OVERWRITTEN``.
    """

    input_file: Optional[str] = None
    output_file: Optional[str] = None
    input_files: tuple[str, ...] = ()
    output_files: tuple[str, ...] = ()
    fingerprints: tuple[FileFingerprint, ...] = ()
    captured_at: Optional[str] = None
    state: ArtifactState = ArtifactState.CURRENT

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
        object.__setattr__(
            self,
            "fingerprints",
            tuple(
                value if isinstance(value, FileFingerprint) else FileFingerprint(**value)
                for value in self.fingerprints
            ),
        )
        if not isinstance(self.state, ArtifactState):
            object.__setattr__(self, "state", ArtifactState(self.state))

    @property
    def manifests(self):
        return _unique((*self.input_files, *self.output_files))


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
    """Persistent identity and last known state of one ASE2SPRKKR run.

    ``state`` describes the complete ASE2SPRKKR lifecycle. ``backend_state``
    records the last status reported by the execution backend. A terminal
    backend state therefore normally corresponds to ``state == FINISHING``
    until parsing, stage-out and local finalization have completed.

    ``finished_at`` records the first transition of the ASE2SPRKKR lifecycle to
    ``FINISHED``, ``FAILED`` or ``CANCELLED`` and is immutable afterwards.
    Later metadata changes, including artifact invalidation, update
    ``updated_at`` but never ``finished_at``.

    Artifact freshness is independent of calculation outcome. A successful run
    remains ``FINISHED`` when its files are later changed; only
    ``artifacts.state`` becomes ``OVERWRITTEN``.
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
    backend_state: Optional[RunState] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    artifacts: Optional[ArtifactRecord] = None
    schema_version: int = _RUN_RECORD_SCHEMA_VERSION

    def __post_init__(self):
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        if isinstance(self.submission, dict):
            object.__setattr__(self, "submission", SubmissionKey(**self.submission))
        if isinstance(self.reference, dict):
            object.__setattr__(self, "reference", RunReference(**self.reference))
        if isinstance(self.artifacts, dict):
            object.__setattr__(self, "artifacts", ArtifactRecord(**self.artifacts))
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
    def _decode(path):
        try:
            with path.open() as file:
                data = json.load(file)
            # ``finished_at`` was added without changing the v1 schema number.
            # Preserve records written by earlier development versions rather
            # than treating an absent optional timestamp as corruption.
            data.setdefault("finished_at", None)
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
        with self._lock, _filesystem_lock(self.directory / ".submit.lock"):
            yield

    def save(self, record: RunRecord):
        """Atomically store a complete record supplied by the caller."""
        with self._lock:
            return self._write_unlocked(record)

    def update(self, record: RunRecord, updater):
        """Read latest, merge under a per-run lock, and atomically replace it."""
        with self._lock, _filesystem_lock(self._run_lock_path(record)):
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
        """Return stored runs without querying backends or touching artifacts."""
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

    @staticmethod
    def _reader_for(record):
        from ..outputs.task_result import KkrOutputReader

        return KkrOutputReader.class_for_task(record.task)(
            None, None, record.directory
        )

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
    def _artifacts_for_spec(spec):
        directory = os.path.abspath(os.fspath(spec.directory))
        input_file = _artifact_path(getattr(spec, "input_file", None), directory)
        output_file = _artifact_path(getattr(spec, "output_file", None), directory)
        input_files = _unique(
            (
                input_file,
                *(
                    _artifact_path(value, directory)
                    for value in getattr(spec, "input_files", ())
                ),
            )
        )
        output_files = _unique(
            (
                output_file,
                *(
                    _artifact_path(value, directory)
                    for value in getattr(spec, "output_files", ())
                ),
            )
        )
        return ArtifactRecord(
            input_file=input_file,
            output_file=output_file,
            input_files=input_files,
            output_files=output_files,
        )

    @staticmethod
    def _artifact_snapshot(artifacts):
        if artifacts is None:
            return ()
        fingerprints = []
        seen = set()
        for manifest in artifacts.manifests:
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
    def _artifacts_missing(artifacts):
        if artifacts is None:
            return False
        for manifest in artifacts.manifests:
            if glob.has_magic(manifest):
                if not glob.glob(manifest, recursive=True):
                    return True
            elif not os.path.isfile(manifest):
                return True
        return False

    def _capture_artifacts(self, artifacts):
        if artifacts is None or artifacts.captured_at is not None:
            return artifacts
        return replace(
            artifacts,
            fingerprints=self._artifact_snapshot(artifacts),
            captured_at=_now(),
            state=ArtifactState.CURRENT,
        )

    def new_record(
        self,
        *,
        task,
        directory,
        dataset=None,
        endpoint=None,
        submission=None,
        artifacts=None,
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
            dataset=None if dataset is None else str(dataset),
            directory=os.path.abspath(os.fspath(directory)),
            endpoint=endpoint,
            submission=submission,
            artifacts=artifacts,
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
    def _merged_artifacts(latest, requested):
        if requested is _UNCHANGED or requested is None:
            return latest.artifacts
        if isinstance(requested, dict):
            requested = ArtifactRecord(**requested)
        current = latest.artifacts
        if current is None:
            return requested
        if current.state == ArtifactState.OVERWRITTEN:
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
        artifacts=_UNCHANGED,
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
            merged_artifacts = self._merged_artifacts(latest, artifacts)
            if merged_artifacts is not latest.artifacts:
                changes["artifacts"] = merged_artifacts
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
        """Persist one backend status without declaring local finalization done."""
        backend_state = RunState(status.state)
        state = RunState.FINISHING if backend_state in _FINAL_STATES else backend_state
        return self.transition(
            record,
            state=state,
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
        """Publish local finalization and capture immutable artifact fingerprints."""
        artifacts = self._capture_artifacts(record.artifacts)
        if backend_state is None:
            backend_state = state if state in _FINAL_STATES else _UNCHANGED
        return self.transition(
            record,
            state=state,
            backend_state=backend_state,
            reference=reference,
            return_code=return_code,
            message=message,
            artifacts=artifacts,
        )

    def check_artifacts(self, run):
        """Refresh only artifact freshness, never the calculation outcome."""
        record = (
            self.store.get(run)
            if isinstance(run, str)
            else self.store.get(run.id)
        )
        artifacts = record.artifacts
        if (
            record.state not in _FINAL_STATES
            or artifacts is None
            or artifacts.captured_at is None
            or artifacts.state == ArtifactState.OVERWRITTEN
        ):
            return record
        if self._artifact_snapshot(artifacts) == artifacts.fingerprints:
            return record

        def mark_overwritten(latest):
            current = latest.artifacts
            if current is None or current.state == ArtifactState.OVERWRITTEN:
                return latest
            return replace(
                latest,
                updated_at=_now(),
                artifacts=replace(current, state=ArtifactState.OVERWRITTEN),
            )

        return self.store.update(record, mark_overwritten)

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
            dataset=dataset,
            directory=spec.directory,
            endpoint=endpoint,
            submission=submission,
            artifacts=self._artifacts_for_spec(spec),
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
        executor=None,
        reader=None,
        callback=None,
        output_listener=None,
    ):
        """Attach to a stored run, reconstructing its persisted endpoint by default.

        ``executor`` is retained as an explicit override for legacy records that
        predate endpoint persistence and for custom backends that cannot be
        reconstructed from :class:`ExecutorEndpoint` alone.

        ``output_listener`` is installed by the backend before monitoring starts.
        Attached runners therefore cannot lose output replayed immediately during
        attach, including output from an already completed run.
        """
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if reader is None:
            reader = self._reader_for(record)
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

    def load_result(self, run, reader=None):
        """Load a historical TaskResult, fetching missing remote artifacts on demand.

        Ordinary history listing remains metadata-only. Network/filesystem
        retrieval is performed here only when a declared artifact is missing
        locally and the persisted endpoint supports ``fetch_artifacts``.
        Missing output after that attempt is represented by an incomplete
        task-specific result. The returned object exposes the refreshed
        persistent record as ``run_record``.
        """
        record = (
            self.store.get(run)
            if isinstance(run, str)
            else self.store.get(run.id)
        )
        if record.state not in _FINAL_STATES:
            raise RuntimeError(f"Run {record.id} is not finished")

        artifacts = record.artifacts
        if (
            self._artifacts_missing(artifacts)
            and record.reference is not None
        ):
            executor = self.executor_for(record)
            fetch_artifacts = getattr(executor, "fetch_artifacts", None)
            if fetch_artifacts is not None:
                fetch_artifacts(record.reference, artifacts)

        record = self.check_artifacts(record)
        artifacts = record.artifacts
        if artifacts is None or artifacts.output_file is None:
            raise RuntimeError(f"Run {record.id} has no persisted output artifact")

        if reader is None:
            reader = self._reader_for(record)

        output_file = artifacts.output_file
        input_file = artifacts.input_file
        return_code = record.return_code
        if return_code is None:
            return_code = 0 if record.state == RunState.FINISHED else 1
        if os.path.isfile(output_file):
            result = reader.read_from_file(
                output_file,
                return_code=return_code,
                input_file=input_file,
                directory=record.directory,
            )
        else:
            result = reader.create_result(output_file, input_file)
            result.complete(None, return_code)
        result.run_record = record
        return result


class PersistentRunner:
    """Caller-owned runner that persists lifecycle transitions."""

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
        state = self._runner.state
        if state == RunState.NEW and self._record.state != RunState.NEW:
            state = self._record.state
        elif state in _FINAL_STATES and self._record.state not in _FINAL_STATES:
            # A backend/job may become terminal before the runner has completed
            # parsing, stage-out and result finalization. Only the completion
            # callback below is allowed to publish the final lifecycle state.
            state = RunState.FINISHING
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

    def submit(self):
        if self._record.state != RunState.NEW:
            raise RuntimeError(f"Run is already {self._record.state.value}")
        self._record = self.manager.reserve(self._record)
        try:
            result = self._runner.submit()
        except BaseException:
            reference = self.reference
            if reference is not None:
                state = self._runner.state
                if state == RunState.NEW:
                    state = RunState.QUEUED
                self._record = self.manager.transition(
                    self._record,
                    state=state,
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
        state = self._runner.state
        if state == RunState.NEW:
            state = RunState.QUEUED
        self._record = self.manager.transition(
            self._record,
            state=state,
            reference=reference,
        )
        return result

    def _completion_state(self):
        backend_state = self._runner.state
        exception = getattr(self._runner, "_exception", None)
        state = (
            RunState.FAILED
            if exception is not None and backend_state != RunState.CANCELLED
            else backend_state
        )
        return state, backend_state, exception

    def _sync_completion(self, result):
        state, backend_state, exception = self._completion_state()
        self._record = self.manager.finalize(
            self._record,
            state=state,
            backend_state=backend_state,
            reference=self.reference,
            return_code=getattr(result, "return_code", None),
            message=None if exception is None else str(exception),
        )

    def wait(self):
        try:
            return self._runner.wait()
        finally:
            if self._record.state != RunState.NEW:
                state, backend_state, exception = self._completion_state()
                self._record = self.manager.finalize(
                    self._record,
                    state=state,
                    backend_state=backend_state,
                    reference=self.reference,
                    return_code=getattr(self.result, "return_code", None),
                    message=None if exception is None else str(exception),
                )

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

    def load_result(self, reader=None):
        """Reload this run from its persisted artifacts."""
        result = self.manager.load_result(self._record, reader=reader)
        self._record = result.run_record
        return result

    def cancel(self):
        try:
            return self._runner.cancel()
        finally:
            if self._record.state != RunState.NEW:
                self._record = self.manager.transition(
                    self._record,
                    state=self._runner.state,
                    reference=self.reference,
                )
