"""Persistent run history and lifecycle orchestration."""

from __future__ import annotations

import asyncio
import json
import os
import warnings
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
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
    """Persistent identity and last known state of one ASE2SPRKKR run."""

    id: str
    state: RunState
    created_at: str
    updated_at: str
    task: str
    directory: str
    dataset: Optional[str] = None
    endpoint: Optional[ExecutorEndpoint] = None
    submission: Optional[SubmissionKey] = None
    reference: Optional[RunReference] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    schema_version: int = _RUN_RECORD_SCHEMA_VERSION

    def __post_init__(self):
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        if isinstance(self.submission, dict):
            object.__setattr__(self, "submission", SubmissionKey(**self.submission))
        if isinstance(self.reference, dict):
            object.__setattr__(self, "reference", RunReference(**self.reference))
        if not isinstance(self.state, RunState):
            object.__setattr__(self, "state", RunState(self.state))
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
                return RunRecord(**json.load(file))
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
        """Return stored runs without querying execution backends."""
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

    def __init__(self, store=None, *, executor_pool=None):
        self.store = store if store is not None else RunStore()
        self.executor_pool = (
            executor_pool if executor_pool is not None else default_executor_pool()
        )
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

    def new_record(
        self,
        *,
        task,
        directory,
        dataset=None,
        endpoint=None,
        submission=None,
        run_id=None,
    ):
        """Create the in-memory identity of an unsubmitted persistent run."""
        now = _now()
        return RunRecord(
            id=run_id or uuid4().hex,
            state=RunState.NEW,
            created_at=now,
            updated_at=now,
            task=str(task).upper(),
            dataset=None if dataset is None else str(dataset),
            directory=os.path.abspath(os.fspath(directory)),
            endpoint=endpoint,
            submission=submission,
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

    def transition(
        self,
        record,
        *,
        state=None,
        reference=_UNCHANGED,
        return_code=_UNCHANGED,
        message=_UNCHANGED,
    ):
        """Merge one lifecycle update into the latest persisted record."""
        def merge(latest):
            final_locked = latest.state in _FINAL_STATES
            changes = {
                "state": self._merged_state(latest, state),
                "updated_at": _now(),
            }
            # Backend identity is monotonic: once known, a stale snapshot must
            # neither erase nor replace it.
            if (
                reference is not _UNCHANGED
                and reference is not None
                and latest.reference is None
            ):
                changes["reference"] = reference

            # Once a final outcome is known, later stale callbacks may only fill
            # missing details. They cannot replace known precise information.
            if return_code is not _UNCHANGED and return_code is not None:
                if not final_locked or latest.return_code is None:
                    changes["return_code"] = return_code
            if message is not _UNCHANGED and message is not None:
                if not final_locked or latest.message is None:
                    changes["message"] = message
            return replace(latest, **changes)

        return self.store.update(record, merge)

    def reserve(self, record):
        """Atomically reserve *record*'s persistent submission namespace."""
        if record.submission is None:
            raise ValueError("Persistent run lacks a submission key")
        with self.store.submission_lock():
            conflicts = tuple(
                active
                for active in self.store.runs(active=True)
                if active.id != record.id
                and active.submission is not None
                and record.submission.conflicts_with(active.submission)
            )
            if conflicts:
                raise DuplicateSubmissionError(record, conflicts)
            return self.transition(
                record,
                state=RunState.SUBMITTING,
                reference=None,
            )

    def update(self, record, status: RunStatus):
        return self.transition(
            record,
            state=status.state,
            reference=status.reference,
            return_code=status.return_code,
            message=status.message,
        )

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

    def attach(self, run, executor, reader, callback=None):
        """Attach to a stored run using an explicitly supplied backend."""
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        endpoint = record.endpoint or endpoint_for_executor(executor)
        self.register_executor(executor, endpoint, record.reference)
        return self._wrap_runner(
            executor,
            record,
            lambda completed: executor.attach(
                record.reference, reader, callback=completed
            ),
            callback,
        )


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

    def _sync_completion(self, result):
        self._record = self.manager.transition(
            self._record,
            state=self._runner.state,
            reference=self.reference,
            return_code=getattr(result, "return_code", None),
        )

    def wait(self):
        try:
            return self._runner.wait()
        finally:
            if self._record.state != RunState.NEW:
                self._record = self.manager.transition(
                    self._record,
                    state=self._runner.state,
                    reference=self.reference,
                    return_code=getattr(self.result, "return_code", None),
                )

    def run(self):
        if self._record.state == RunState.NEW:
            self.submit()
        return self.wait()

    async def run_async(self):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.run)

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
