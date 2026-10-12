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


_RUN_RECORD_SCHEMA_VERSION = 1
_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}
_UNCHANGED = object()


def _now():
    return datetime.now(timezone.utc).isoformat()


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def _endpoint_identity(executor):
    """Return a stable, serializable identity of the backend connection."""
    descriptor = {
        "executor": executor.name
        or f"{type(executor).__module__}.{type(executor).__qualname__}"
    }
    for name in ("backend", "url", "computer", "command", "backend_options"):
        if hasattr(executor, name):
            descriptor[name] = getattr(executor, name)
    transport = getattr(executor, "transport", None)
    if transport is not None:
        transport_descriptor = {
            "type": f"{type(transport).__module__}.{type(transport).__qualname__}"
        }
        for name in ("host", "ssh", "options"):
            if hasattr(transport, name):
                transport_descriptor[name] = getattr(transport, name)
        descriptor["transport"] = transport_descriptor
    return json.dumps(_jsonable(descriptor), sort_keys=True, separators=(",", ":"))


class CorruptRunRecord(RuntimeError):
    """A persisted run record cannot be decoded or validated."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Corrupt run record {self.path}: {message}")


@dataclass(frozen=True)
class SubmissionKey:
    """Persistent namespace used only for duplicate-submit protection."""

    endpoint: str
    directory: str
    task: str
    dataset: Optional[str] = None

    def __post_init__(self):
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
    submission: Optional[SubmissionKey] = None
    reference: Optional[RunReference] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    schema_version: int = _RUN_RECORD_SCHEMA_VERSION

    def __post_init__(self):
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

    @contextmanager
    def submission_lock(self):
        """Serialize duplicate checking and the SUBMITTING record creation."""
        self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        path = self.directory / ".submit.lock"
        with self._lock, path.open("a+b") as lock_file:
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

    def save(self, record: RunRecord):
        with self._lock:
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
    :class:`PersistentRunner` objects remain owned by their callers.
    """

    def __init__(self, store=None):
        self.store = store if store is not None else RunStore()
        self._executors = {}
        self._lock = RLock()

    @staticmethod
    def _executor_key(reference):
        return reference.executor, reference.computer

    def register_executor(self, executor: Executor, reference=None):
        """Remember a live backend only for status refresh in this process."""
        if reference is not None:
            key = self._executor_key(reference)
        else:
            name = executor.name
            computer = getattr(executor, "computer", None)
            if name == "psij" and hasattr(executor, "backend"):
                name = f"psij:{executor.backend}"
                computer = getattr(executor, "url", None)
            elif name == "remote" and hasattr(executor, "backend"):
                name = f"remote:{executor.backend}"
            key = name, computer
        with self._lock:
            self._executors[key] = executor
        return executor

    def executor_for(self, reference):
        with self._lock:
            executor = self._executors.get(self._executor_key(reference))
            if executor is None and reference.computer is not None:
                executor = self._executors.get((reference.executor, None))
        if executor is None:
            raise KeyError(
                f"No live executor registered for {reference.executor} on "
                f"{reference.computer or 'local machine'}"
            )
        return executor

    def new_record(
        self, *, task, directory, dataset=None, submission=None, run_id=None
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
            submission=submission,
        )

    def transition(
        self,
        record,
        *,
        state=None,
        reference=_UNCHANGED,
        return_code=_UNCHANGED,
        message=_UNCHANGED,
    ):
        """Persist one monotonic information update for *record*."""
        new_state = record.state if state is None else RunState(state)
        if record.state in _FINAL_STATES and new_state not in _FINAL_STATES:
            new_state = record.state

        changes = {"state": new_state, "updated_at": _now()}
        if reference is not _UNCHANGED:
            changes["reference"] = reference
        if return_code is not _UNCHANGED and return_code is not None:
            changes["return_code"] = return_code
        if message is not _UNCHANGED and message is not None:
            changes["message"] = message
        return self.store.save(replace(record, **changes))

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
        submission = SubmissionKey(
            endpoint=_endpoint_identity(executor),
            directory=executor.submission_directory(spec),
            task=task,
            dataset=dataset,
        )
        record = self.new_record(
            task=task,
            dataset=dataset,
            directory=spec.directory,
            submission=submission,
            run_id=run_id,
        )
        self.register_executor(executor)
        return self._wrap_runner(
            executor,
            record,
            lambda completed: executor.create_run(reader, spec, callback=completed),
            callback,
        )

    def refresh(self, records=None):
        """Refresh referenced runs, batching each live executor into one query."""
        records = tuple(self.store.runs(active=True) if records is None else records)
        current = {record.id: record for record in records}
        groups = {}
        for record in records:
            if record.reference is None:
                continue
            executor = self.executor_for(record.reference)
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
        self.register_executor(executor, record.reference)
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
        self.manager.register_executor(self.executor, reference)
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
