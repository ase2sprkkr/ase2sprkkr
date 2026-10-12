from dataclasses import asdict
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, Executor, RunReference, RunState, RunTarget
from ..remote import RemoteExecutor
from ..tracking import (
    RunConflictError,
    RunRegistry,
    TrackingExecutor,
    execution_target,
)


class _Runner:
    def __init__(self, reference=None):
        self.result = object()
        self._state = RunState.NEW
        self._reference = reference
        self.submit_count = 0

    @property
    def state(self):
        return self._state

    @property
    def reference(self):
        return self._reference

    def submit(self):
        if self._state != RunState.NEW:
            raise RuntimeError("already submitted")
        self.submit_count += 1
        self._state = RunState.RUNNING

    def wait(self):
        self._state = RunState.FINISHED
        return self.result

    def cancel(self):
        self._state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass


class _Executor(Executor):
    name = "test"

    def __init__(self):
        self.runners = []

    def create_run(self, reader, spec, callback=None):
        runner = _Runner()
        self.runners.append(runner)
        return runner


def _spec(tmp_path, *, input_name="input.inp", output_name="output.out", cache=()):
    return ExecutionSpec(
        argv=("true",),
        input_file=str(tmp_path / input_name),
        output_file=str(tmp_path / output_name),
        directory=str(tmp_path),
        input_files=(str(tmp_path / input_name),),
        output_files=(str(tmp_path / output_name),),
        cache_files=tuple(str(tmp_path / name) for name in cache),
    )


def test_run_target_conflicts_only_on_shared_writes():
    first = RunTarget("local", ("shared.in",), ("first.out",))
    same_input = RunTarget("local", ("shared.in",), ("second.out",))
    reads_first_output = RunTarget("local", ("first.out",), ("third.out",))
    other_machine = RunTarget("ssh:other", (), ("first.out",))

    assert not first.conflicts_with(same_input)
    assert first.conflicts_with(reads_first_output)
    assert not first.conflicts_with(other_machine)


def test_execution_target_treats_cache_as_read_write(tmp_path):
    target = Executor().execution_target(_spec(tmp_path, cache=("cache.dat",)))
    cache = str((tmp_path / "cache.dat").resolve())

    assert cache in target.read_files
    assert cache in target.write_files


def test_registry_reserves_atomically_and_releases_on_wait(tmp_path):
    registry = RunRegistry()
    executor = TrackingExecutor(_Executor(), registry)
    first = executor.create_run(None, _spec(tmp_path))
    second = executor.create_run(None, _spec(tmp_path))

    first.submit()
    assert second.conflicts == (first,)
    with pytest.raises(RunConflictError) as exc:
        second.submit()
    assert exc.value.conflicts == (first,)

    first.wait()
    assert second.conflicts == ()
    second.submit()


def test_registry_allows_explicit_conflict_override(tmp_path):
    registry = RunRegistry()
    executor = TrackingExecutor(_Executor(), registry)
    first = executor.create_run(None, _spec(tmp_path))
    second = executor.create_run(None, _spec(tmp_path))

    first.submit()
    second.submit(allow_conflicts=True)

    assert set(registry.active) == {first, second}


def test_different_outputs_in_same_directory_can_run_together(tmp_path):
    registry = RunRegistry()
    executor = TrackingExecutor(_Executor(), registry)
    first = executor.create_run(None, _spec(tmp_path, output_name="one.out"))
    second = executor.create_run(None, _spec(tmp_path, output_name="two.out"))

    first.submit()
    second.submit()

    assert set(registry.active) == {first, second}


def test_run_reference_carries_serializable_target():
    target = RunTarget("ssh:cluster", ("/work/in",), ("/work/out",))
    reference = RunReference("remote:local", "42", target=target)
    data = asdict(reference)

    restored = RunReference(**data)
    assert restored.target == target


class _RemoteTransport:
    host = "cluster-a"

    async def run(self, argv, **kwargs):
        return SimpleNamespace(stdout=b"/home/test")


def test_remote_target_uses_remote_paths_and_machine_identity(tmp_path):
    executor = RemoteExecutor(
        _RemoteTransport(),
        local_root=tmp_path,
        remote_root="/scratch/work",
    )
    spec = _spec(tmp_path)

    target = execution_target(executor, spec)

    assert target.location == "remote:cluster-a"
    assert "/scratch/work/input.inp" in target.read_files
    assert "/scratch/work/output.out" in target.write_files
    assert str(tmp_path / "output.out") not in target.write_files
