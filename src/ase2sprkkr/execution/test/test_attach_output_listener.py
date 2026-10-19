from dataclasses import replace
from types import SimpleNamespace

from ase2sprkkr.execution.base import LaunchSpec, RunState
from ase2sprkkr.execution.manager import RunRecord
from ase2sprkkr.execution.run_model import Execution, WorkingDirectory, FileManifest, FileEntry, FileMapping
from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.local import LocalExecutor
from ase2sprkkr.execution.manager import ExecutionManager, RunStore
from ase2sprkkr.execution.psij import PsijExecutor, PsijRunner
from ase2sprkkr.execution.remote_endpoint import _RemoteRun, _RemoteRunContext


class _AttachedRunner:
    def __init__(self, run):
        self.backend_id = run.execution.backend_id
        self.result = SimpleNamespace(return_code=0)
        self.state = RunState.RUNNING
        self.exception = None
        self.failure_phase = None
        self.failure_message = None

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass

    def start(self):
        self.state = RunState.FINISHED

    def wait(self):
        return self.result

    def cancel(self):
        pass


class _AttachExecutor(LocalExecutor):
    def acknowledge(self, run):
        # Synthetic backend IDs in this focused listener test are not real
        # supervisor handles; acknowledgment is tested separately.
        pass

    def attach(self, run, callback=None, output_listener=None, start=True):
        if output_listener is not None:
            output_listener("manager replay\n", "out")
        runner = _AttachedRunner(run)
        if start:
            runner.start()
        return runner


def test_manager_attach_forwards_initial_output_listener(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(task="SCF", directory=tmp_path)
    record = replace(
        record, state=RunState.RUNNING,
        execution=record.execution.submitted("123"),
    )
    manager.store.save(record)

    output = []
    runner = manager.attach(
        record,
        executor=_AttachExecutor(),
        process_result=lambda result: result,
        output_listener=lambda text, kind: output.append((text, kind)),
    )
    runner.wait()

    assert output == [("manager replay\n", "out")]


def _remote_context(tmp_path):
    class Endpoint:
        backend = "local"
        backend_options = {}
        command = ("ase2sprkkr", "remote")
        computer = "host"
        transport = SimpleNamespace()

        async def _remote_home(self):
            return "/remote"

        async def _status(self, directory, run_id, files=None):
            return {
                "state": "running",
                "input_file": "/remote/work/input.inp",
                "output_file": "/remote/work/output.out",
            }

    context = _RemoteRunContext(
        Endpoint(),
        LaunchSpec(local_root=str(tmp_path), remote_root="/remote"),
    )
    run = RunRecord(
        id="run", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(executor="ssh", host="host"),
            WorkingDirectory(str(tmp_path / "work"), "/remote/work"), "run",
        ),
        files=FileManifest(
            stdout=FileEntry(
                local_path=str(tmp_path / "work" / "output.out"),
                remote_path="/remote/work/output.out",
            ),
            mapping=FileMapping(str(tmp_path), "/remote"),
        ),
    )
    return context, run


def test_remote_attach_installs_listener_before_monitor(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"remote replay\n", "out")
        if self.callback is not None:
            self.callback(self.result)
        self._done.set()

    monkeypatch.setattr(_RemoteRun, "_thread_main", immediate_monitor)
    context, run = _remote_context(tmp_path)

    import asyncio

    runner = asyncio.run(
        context.attach(
            run,
            callback=lambda result: events.append(("finished", None)),
            output_listener=lambda text, kind: events.append((text, kind)),
        )
    )
    runner._done.wait(1)

    assert events == [("remote replay\n", "out"), ("finished", None)]


def test_remote_attach_can_defer_monitor_start(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"remote replay\n", "out")
        self._done.set()

    monkeypatch.setattr(_RemoteRun, "_thread_main", immediate_monitor)
    context, run = _remote_context(tmp_path)

    import asyncio

    runner = asyncio.run(
        context.attach(
            run,
            output_listener=lambda text, kind: events.append((text, kind)),
            start=False,
        )
    )

    assert runner._thread is None
    assert events == []

    runner.start()
    runner._done.wait(1)
    assert events == [("remote replay\n", "out")]


def _psij_executor():
    class Job:
        def __init__(self):
            self.native_id = "42"
            self.status = SimpleNamespace(state=SimpleNamespace(name="ACTIVE"))

    class JobExecutor:
        def attach(self, job, native_id):
            job.native_id = native_id

    executor = object.__new__(PsijExecutor)
    executor.backend = "slurm"
    executor.url = None
    executor.poll_interval = 0.01
    executor._psij = SimpleNamespace(Job=Job)
    executor._job_executor = JobExecutor()
    executor._status_jobs = {}
    return executor


def _psij_run(tmp_path):
    return RunRecord(
        id="psij-run", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(executor="psij", psij_backend="slurm"),
            WorkingDirectory(str(tmp_path)), "42",
        ),
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "output.out")),
            stderr=FileEntry(local_path=str(tmp_path / "error.err")),
        ),
    )


def test_psij_attach_installs_listener_before_monitor(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"psij replay\n", "out")
        if self.callback is not None:
            self.callback(self.result)
        self._done.set()

    monkeypatch.setattr(PsijRunner, "_thread_main", immediate_monitor)

    runner = _psij_executor().attach(
        _psij_run(tmp_path),
        callback=lambda result: events.append(("finished", None)),
        output_listener=lambda text, kind: events.append((text, kind)),
    )
    runner._done.wait(1)

    assert events == [("psij replay\n", "out"), ("finished", None)]


def test_psij_attach_can_defer_monitor_start(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"psij replay\n", "out")
        self._done.set()

    monkeypatch.setattr(PsijRunner, "_thread_main", immediate_monitor)

    runner = _psij_executor().attach(
        _psij_run(tmp_path),
        output_listener=lambda text, kind: events.append((text, kind)),
        start=False,
    )

    assert runner._thread is None
    assert events == []

    runner.start()
    runner._done.wait(1)
    assert events == [("psij replay\n", "out")]
