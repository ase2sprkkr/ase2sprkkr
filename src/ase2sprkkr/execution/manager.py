"""Persistent bookkeeping and bulk refresh for managed execution runs."""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from threading import RLock
from typing import Optional
from uuid import uuid4

from platformdirs import user_state_dir

from .base import Executor, RunReference, RunState, RunStatus


_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED, RunState.LOST}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _identity_backend(executor):
    """Unwrap executor decorators only for persistent backend identity."""
    backend = executor
    seen = set()
    while id(backend) not in seen:
        seen.add(id(backend))
        name = getattr(backend, "name", None)
        if name not in ("psij", "remote") or hasattr(backend, "backend"):
            break
        nested = getattr(backend, "executor", None)
        if nested is None:
            break
        backend = nested
    return backend


@dataclass(frozen=True)
class RunRecord:
    """Persistent ase2sprkkr identity and last known state of one run."""

    id: str
    reference: RunReference
    state: RunState
    created_at: str
    updated_at: str
    return_code: Optional[int] = None
    message: Optional[str] = None

    def __post_init__(self):
        if isinstance(self.reference, dict):
            object.__setattr__(self, "reference", RunReference(**self.reference))
        if not isinstance(self.state, RunState):
            object.__setattr__(self, "state", RunState(self.state))


class RunStore:
    """Small persistent store using one atomic JSON file per run."""

    def __init__(self, directory=None):
        if directory is None:
            directory = Path(user_state_dir("ase2sprkkr")) / "runs"
        self.directory = Path(directory).expanduser()
        self._lock = RLock()

    def _path(self, run_id):
        return self.directory / (str(run_id) + ".json")

    def save(self, record: RunRecord):
        with self._lock:
            self.directory.mkdir(parents=True, exist_ok=True)
            path = self._path(record.id)
            temporary = path.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(asdict(record), sort_keys=True))
            os.replace(temporary, path)
        return record

    def get(self, run_id):
        path = self._path(run_id)
        with path.open() as file:
            return RunRecord(**json.load(file))

    def runs(self, *, active=None):
        with self._lock:
            if not self.directory.is_dir():
                return ()
            records = []
            for path in self.directory.glob("*.json"):
                with path.open() as file:
                    record = RunRecord(**json.load(file))
                if active is True and record.state in _FINAL_STATES:
                    continue
                if active is False and record.state not in _FINAL_STATES:
                    continue
                records.append(record)
            return tuple(sorted(records, key=lambda record: record.created_at))


class ExecutorPool:
    """Explicit application-owned cache for long-lived executor instances."""

    def __init__(self):
        self._executors = {}
        self._lock = RLock()

    def get(self, key, factory):
        with self._lock:
            executor = self._executors.get(key)
            if executor is None:
                executor = factory()
                self._executors[key] = executor
            return executor

    def for_profile(self, profile):
        """Return one shared backend for an immutable SPR-KKR runtime profile."""
        from ..sprkkr.profile import executor_from_profile

        return self.get(profile, lambda: executor_from_profile(profile))


class RunManager:
    """Connect persistent ase2sprkkr run records with execution backends."""

    def __init__(self, store=None, executor_pool=None):
        self.store = store if store is not None else RunStore()
        self.executor_pool = executor_pool if executor_pool is not None else ExecutorPool()
        self._executors = {}
        self._managed_executors = {}
        self._lock = RLock()

    @staticmethod
    def _executor_key(reference):
        return reference.executor, reference.computer

    def register_executor(self, executor: Executor, *, executor_name=None, computer=None):
        """Register a live backend or decorator for persistent references."""
        identity = _identity_backend(executor)
        name = executor_name or getattr(identity, "name", executor.name)
        if name == "psij" and hasattr(identity, "backend"):
            name = "psij:{}".format(identity.backend)
            if computer is None:
                computer = getattr(identity, "url", None)
        elif name == "remote" and hasattr(identity, "backend"):
            name = "remote:{}".format(identity.backend)
            if computer is None:
                computer = getattr(identity, "computer", None)
        with self._lock:
            self._executors[(name, computer)] = executor
        return executor

    def executor_for(self, reference):
        """Return the registered backend capable of refreshing *reference*."""
        with self._lock:
            executor = self._executors.get(self._executor_key(reference))
            if executor is None and reference.computer is not None:
                executor = self._executors.get((reference.executor, None))
        if executor is None:
            raise KeyError(
                "No executor registered for {} on {}".format(
                    reference.executor, reference.computer or "local machine"
                )
            )
        return executor

    def backend_for_profile(self, profile):
        """Return and register the shared bare backend for *profile*."""
        return self.register_executor(self.executor_pool.for_profile(profile))

    def executor_for_profile(self, profile):
        """Return one managed executor over the shared backend for *profile*."""
        with self._lock:
            managed = self._managed_executors.get(profile)
            if managed is not None:
                return managed
        backend = self.backend_for_profile(profile)
        managed = ManagedExecutor(backend, self)
        with self._lock:
            existing = self._managed_executors.setdefault(profile, managed)
        return existing

    def remember(self, reference, *, state=RunState.QUEUED, run_id=None):
        if reference is None:
            return None
        now = _now()
        record = RunRecord(
            id=run_id or uuid4().hex,
            reference=reference,
            state=RunState(state),
            created_at=now,
            updated_at=now,
        )
        return self.store.save(record)

    def update(self, record, status: RunStatus):
        updated = replace(
            record,
            state=status.state,
            updated_at=_now(),
            return_code=status.return_code,
            message=status.message,
        )
        return self.store.save(updated)

    def refresh(self, records=None):
        """Refresh active records, batching each executor identity into one query."""
        if records is None:
            records = self.store.runs(active=True)
        groups = {}
        for record in records:
            executor = self.executor_for(record.reference)
            groups.setdefault(executor, []).append(record)

        updated = []
        for executor, group in groups.items():
            statuses = executor.status(*(record.reference for record in group))
            if len(statuses) != len(group):
                raise RuntimeError("Executor returned wrong number of run statuses")
            updated.extend(
                self.update(record, status)
                for record, status in zip(group, statuses)
            )
        return tuple(updated)

    def attach(self, run, reader, callback=None):
        """Reconnect to a stored persistent run and keep updating its record."""
        record = self.store.get(run) if isinstance(run, str) else run
        executor = self.executor_for(record.reference)
        runner = executor.attach(record.reference, reader, callback=callback)
        return ManagedRunner(self, runner, record=record)


class ManagedRunner:
    """Runner proxy that persists its reference and last known state."""

    def __init__(self, manager, runner, *, record=None):
        self.manager = manager
        self._runner = runner
        self._record = record

    @property
    def result(self):
        return self._runner.result

    @property
    def reference(self):
        return self._runner.reference

    @property
    def state(self):
        state = self._runner.state
        if self._record is not None and state != self._record.state:
            self._record = self.manager.update(
                self._record, RunStatus(self.reference, state)
            )
        return state

    def add_output_listener(self, callback):
        return self._runner.add_output_listener(callback)

    def submit(self, *args, **kwargs):
        out = self._runner.submit(*args, **kwargs)
        reference = self._runner.reference
        if reference is not None and self._record is None:
            self._record = self.manager.remember(reference, state=self._runner.state)
        return out

    def wait(self):
        try:
            return self._runner.wait()
        finally:
            if self._record is not None:
                self._record = self.manager.update(
                    self._record, RunStatus(self.reference, self._runner.state)
                )

    def run(self):
        if self._runner.state == RunState.NEW:
            self.submit()
        return self.wait()

    def cancel(self):
        try:
            return self._runner.cancel()
        finally:
            if self._record is not None:
                self._record = self.manager.update(
                    self._record, RunStatus(self.reference, self._runner.state)
                )


class ManagedExecutor(Executor):
    """Executor decorator that persists submitted run references."""

    def __init__(self, executor, manager: RunManager):
        self.executor = executor
        self.manager = manager
        self.manager.register_executor(executor)

    @property
    def name(self):
        return self.executor.name

    def execution_target(self, spec):
        return self.executor.execution_target(spec)

    def resolve_mpi(self, mpi=None):
        return self.executor.resolve_mpi(mpi)

    def status(self, *references):
        return self.executor.status(*references)

    def create_run(self, reader, spec, callback=None):
        return ManagedRunner(
            self.manager,
            self.executor.create_run(reader, spec, callback=callback),
        )

    def attach(self, reference, reader, callback=None):
        return ManagedRunner(
            self.manager,
            self.executor.attach(reference, reader, callback=callback),
        )
