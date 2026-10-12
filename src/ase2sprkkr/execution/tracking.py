"""Execution-level tracking and collision prevention for active runs."""

from __future__ import annotations

import asyncio
from dataclasses import replace
from functools import singledispatch
from threading import Lock, RLock, Thread

from .base import Executor, RunState, RunTarget
from .remote import RemoteExecutor
from .transport import LocalTransport, SshTransport


_FINAL_STATES = {RunState.FINISHED, RunState.FAILED, RunState.CANCELLED}


def _run_sync(coro):
    """Run a small target-resolution coroutine from synchronous API code."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)

    result = []
    error = []

    def worker():
        try:
            result.append(asyncio.run(coro))
        except BaseException as exc:
            error.append(exc)

    thread = Thread(target=worker)
    thread.start()
    thread.join()
    if error:
        raise error[0]
    return result[0]


@singledispatch
def execution_target(executor: Executor, spec) -> RunTarget:
    """Resolve the concrete filesystem resources owned by ``executor``."""
    return executor.execution_target(spec)


@execution_target.register
def _(executor: RemoteExecutor, spec) -> RunTarget:
    local = Executor.execution_target(executor, spec)

    async def mapped():
        return tuple([await executor.map_path(path) for path in local.files])

    files = _run_sync(mapped())
    transport = executor.transport
    if isinstance(transport, LocalTransport):
        location = "local"
    elif isinstance(transport, SshTransport):
        location = "ssh:{}".format(transport.host)
    else:
        computer = executor.computer
        if computer is None:
            raise ValueError(
                "Remote transport must identify its target machine for run tracking"
            )
        location = "remote:{}".format(computer)
    return RunTarget(location, files)


class RunConflictError(RuntimeError):
    """Raised when a run would own resources already used by an active run."""

    def __init__(self, target, conflicts):
        self.target = target
        self.conflicts = tuple(conflicts)
        super().__init__(
            "Execution target conflicts with {} active run{}".format(
                len(self.conflicts), "" if len(self.conflicts) == 1 else "s"
            )
        )


class TrackedRunner:
    """Runner proxy whose submission is reserved in a :class:`RunRegistry`."""

    def __init__(self, registry, runner, target: RunTarget, *, submitted=False):
        self._registry = registry
        self._runner = runner
        self.target = target
        self._submit_lock = Lock()
        self._submitted = submitted

    @property
    def result(self):
        return self._runner.result

    @property
    def state(self):
        state = self._runner.state
        if state in _FINAL_STATES:
            self._registry._release(self)
        return state

    @property
    def reference(self):
        reference = self._runner.reference
        if reference is None or reference.target == self.target:
            return reference
        return replace(reference, target=self.target)

    @property
    def conflicts(self):
        """Active runs that currently conflict with this run target."""
        return self._registry.conflicts(self.target)

    def add_output_listener(self, callback):
        return self._runner.add_output_listener(callback)

    def submit(self, *, allow_conflicts=False):
        """Submit after atomically reserving this run's resources.

        ``allow_conflicts=True`` is an explicit escape hatch for callers that
        deliberately choose to run despite a previously reported conflict.
        """
        with self._submit_lock:
            if self._submitted:
                return self._runner.submit()

            self._registry._reserve(self, allow_conflicts=allow_conflicts)
            try:
                self._runner.submit()
            except BaseException:
                self._registry._release(self)
                raise
            self._submitted = True

    def _release_if_final(self):
        try:
            state = self._runner.state
        except BaseException:
            return
        if state in _FINAL_STATES:
            self._registry._release(self)

    def wait(self):
        try:
            return self._runner.wait()
        finally:
            if self._submitted:
                self._release_if_final()

    def run(self):
        if not self._submitted:
            self.submit()
        return self.wait()

    async def run_async(self):
        if not self._submitted:
            self.submit()
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.wait)

    def cancel(self):
        try:
            return self._runner.cancel()
        finally:
            self._release_if_final()


class RunRegistry:
    """Thread-safe registry of submitted runs and their reserved resources.

    Registration happens immediately before the backend ``submit()`` call, so
    two callers cannot both pass collision checking before either is visible.
    Failed submissions release their reservation automatically.
    """

    def __init__(self):
        self._lock = RLock()
        self._runs = set()

    def track(self, runner, target: RunTarget):
        """Wrap an unsubmitted runner in this registry."""
        return TrackedRunner(self, runner, target)

    def adopt(self, runner, target: RunTarget):
        """Track a runner that was already submitted before it was attached."""
        run = TrackedRunner(self, runner, target, submitted=True)
        self._reserve(run, allow_conflicts=True)
        return run

    def _prune_locked(self):
        for run in tuple(self._runs):
            try:
                state = run._runner.state
            except BaseException:
                continue
            if state in _FINAL_STATES:
                self._runs.discard(run)

    def conflicts(self, target: RunTarget):
        """Return active tracked runs conflicting with ``target``."""
        with self._lock:
            self._prune_locked()
            return tuple(
                run for run in self._runs if run.target.conflicts_with(target)
            )

    @property
    def active(self):
        """Return currently active tracked runs."""
        with self._lock:
            self._prune_locked()
            return tuple(self._runs)

    def _reserve(self, run, *, allow_conflicts=False):
        with self._lock:
            self._prune_locked()
            conflicts = tuple(
                active
                for active in self._runs
                if active is not run and active.target.conflicts_with(run.target)
            )
            if conflicts and not allow_conflicts:
                raise RunConflictError(run.target, conflicts)
            self._runs.add(run)

    def _release(self, run):
        with self._lock:
            self._runs.discard(run)


class TrackingExecutor(Executor):
    """Executor decorator adding :class:`RunRegistry` tracking to a backend."""

    def __init__(self, executor: Executor, registry: RunRegistry):
        self.executor = executor
        self.registry = registry

    @property
    def name(self):
        return self.executor.name

    def execution_target(self, spec):
        return execution_target(self.executor, spec)

    def resolve_mpi(self, mpi=None):
        return self.executor.resolve_mpi(mpi)

    def create_run(self, reader, spec, callback=None):
        target = self.execution_target(spec)
        runner = self.executor.create_run(reader, spec, callback=callback)
        return self.registry.track(runner, target)

    def attach(self, reference, reader, callback=None):
        runner = self.executor.attach(reference, reader, callback=callback)
        if reference.target is None:
            return runner
        return self.registry.adopt(runner, reference.target)
