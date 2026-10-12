import argparse
import asyncio
import json
import sys
from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionSpec,
    Executor,
    ExecutorPool,
    LocalExecutor,
    ManagedExecutor,
    RunManager,
    RunReference,
    RunState,
    RunStatus,
    RunStore,
    RunTarget,
)
from ase2sprkkr.execution import psij as psij_backend
from ase2sprkkr.execution.remote import _run_directory, _write_state
from ase2sprkkr.tools.commands import remote as remote_command


class _StatusExecutor(Executor):
    name = "test"

    def __init__(self):
        self.calls = []

    def status(self, *references):
        self.calls.append(references)
        return tuple(
            RunStatus(reference, RunState.RUNNING) for reference in references
        )


def test_run_store_roundtrip_and_manager_batches_refresh(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = RunManager(store=store)
    executor = _StatusExecutor()
    manager.register_executor(executor)

    first = manager.remember(
        RunReference(
            executor="test",
            native_id="one",
            target=RunTarget("local", ("one.out",)),
        )
    )
    second = manager.remember(
        RunReference(executor="test", native_id="two")
    )

    restored = RunStore(tmp_path / "runs")
    assert restored.get(first.id).reference.target == RunTarget(
        "local", ("one.out",)
    )

    updated = manager.refresh()

    assert len(executor.calls) == 1
    assert executor.calls[0] == (first.reference, second.reference)
    assert [record.state for record in updated] == [
        RunState.RUNNING,
        RunState.RUNNING,
    ]


def test_executor_pool_reuses_profile_executor():
    pool = ExecutorPool()
    created = []

    def factory():
        value = object()
        created.append(value)
        return value

    first = pool.get(("psij", "slurm"), factory)
    second = pool.get(("psij", "slurm"), factory)

    assert first is second
    assert created == [first]


class _LocalResult:
    pass


class _LocalReader:
    print_output = False

    def create_result(self, output_file, input_file):
        return _LocalResult()

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.return_code = return_code
        return result


def test_managed_local_run_persists_pid_before_submit_returns(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = RunManager(store=store)
    executor = ManagedExecutor(LocalExecutor(), manager)
    output = open(tmp_path / "run.out", "w+b")
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", "import time; time.sleep(1)"),
        input_file=None,
        output_file=output,
        directory=str(tmp_path),
    )
    runner = executor.create_run(_LocalReader(), spec)

    runner.submit()
    reference = runner.reference

    assert reference is not None
    assert reference.executor == "local"
    assert reference.native_id
    records = store.runs(active=True)
    assert len(records) == 1
    assert records[0].reference == reference
    assert LocalExecutor().status(reference)[0].state == RunState.RUNNING

    runner.cancel()
    try:
        runner.wait()
    except RuntimeError:
        pass


class _PsijStatus:
    def __init__(self, name="NEW", exit_code=None, message=None):
        self.state = SimpleNamespace(name=name)
        self.exit_code = exit_code
        self.message = message


class _PsijJob:
    def __init__(self, spec=None):
        self.spec = spec
        self.native_id = None
        self.status = _PsijStatus()


class _PsijBackend:
    def __init__(self):
        self.attached = []

    def attach(self, job, native_id):
        self.attached.append(native_id)
        job.native_id = native_id
        job.status = _PsijStatus("ACTIVE")


class _PsijJobExecutor:
    @staticmethod
    def get_instance(*args, **kwargs):
        raise AssertionError("A supplied backend should be used")


class _PsijJobSpec:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _PsijResourceSpec:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def test_psij_status_reuses_attached_jobs(monkeypatch):
    fake_psij = SimpleNamespace(
        Job=_PsijJob,
        JobSpec=_PsijJobSpec,
        JobExecutor=_PsijJobExecutor,
        ResourceSpecV1=_PsijResourceSpec,
    )
    monkeypatch.setattr(psij_backend, "_load_psij", lambda: fake_psij)
    backend = _PsijBackend()
    executor = psij_backend.PsijExecutor("slurm", job_executor=backend)
    references = (
        RunReference(executor="psij:slurm", native_id="101"),
        RunReference(executor="psij:slurm", native_id="102"),
    )

    first = executor.status(*references)
    second = executor.status(*references)

    assert [status.state for status in first] == [RunState.RUNNING, RunState.RUNNING]
    assert [status.state for status in second] == [RunState.RUNNING, RunState.RUNNING]
    assert backend.attached == ["101", "102"]


def test_remote_status_command_accepts_multiple_run_ids(tmp_path, capsys):
    workspace = tmp_path / "workspace"
    for run_id, state in (("one", "running"), ("two", "finished")):
        _run_directory(workspace, run_id).mkdir(parents=True)
        _write_state(
            workspace,
            run_id,
            id=run_id,
            state=state,
            directory=str(workspace),
        )

    parser = argparse.ArgumentParser()
    remote_command.parser(parser)
    args = parser.parse_args(["status", str(workspace), "one", "two"])
    remote_command.run(args, None)

    states = json.loads(capsys.readouterr().out)
    assert [state["id"] for state in states] == ["one", "two"]
    assert [state["state"] for state in states] == ["running", "finished"]
