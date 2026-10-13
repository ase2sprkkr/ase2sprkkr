from dataclasses import replace
from types import SimpleNamespace

from ase2sprkkr.execution.base import ExecutionSpec, LaunchSpec, RunReference, RunState
from ase2sprkkr.execution.local import LocalExecutor
from ase2sprkkr.execution.manager import ExecutionManager, RunStore
from ase2sprkkr.execution.psij import PsijExecutor, PsijRunner
from ase2sprkkr.execution.remote_endpoint import _RemoteRun, _RemoteRunContext


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None)


class _AttachedRunner:
    def __init__(self, reference, callback=None):
        self.reference = reference
        self.callback = callback
        self.result = SimpleNamespace(return_code=None)
        self.state = RunState.RUNNING

    def add_output_listener(self, callback):
        pass

    def wait(self):
        return self.result

    def cancel(self):
        pass


class _AttachExecutor(LocalExecutor):
    def attach(self, reference, reader, callback=None, output_listener=None):
        if output_listener is not None:
            output_listener("manager replay\n", "out")
        return _AttachedRunner(reference, callback=callback)


def test_manager_attach_forwards_initial_output_listener(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    reference = RunReference(
        executor="local",
        native_id="123",
        directory=str(tmp_path),
    )
    record = manager.new_record(task="SCF", directory=tmp_path)
    record = replace(record, state=RunState.RUNNING, reference=reference)
    manager.store.save(record)

    output = []
    manager.attach(
        record,
        executor=_AttachExecutor(),
        reader=_Reader(),
        output_listener=lambda text, kind: output.append((text, kind)),
    )

    assert output == [("manager replay\n", "out")]


def test_remote_attach_installs_listener_before_monitor(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"remote replay\n", "out")
        if self.callback is not None:
            self.callback(self.result)
        self._done.set()

    monkeypatch.setattr(_RemoteRun, "_thread_main", immediate_monitor)

    class Endpoint:
        backend = "local"
        backend_options = {}
        command = ("ase2sprkkr", "remote")
        computer = "host"
        transport = SimpleNamespace()

        async def _remote_home(self):
            return "/remote"

        async def _status(self, directory, run_id):
            return {
                "state": "running",
                "input_file": "/remote/work/input.inp",
                "output_file": "/remote/work/output.out",
            }

    context = _RemoteRunContext(
        Endpoint(),
        LaunchSpec(local_root=str(tmp_path), remote_root="/remote"),
    )
    reference = RunReference(
        executor="remote:local",
        native_id="run",
        computer="host",
        directory="/remote/work",
        input_file="/remote/work/input.inp",
        output_file="/remote/work/output.out",
    )
    request = {
        "executor": "local",
        "client_sync": {
            "output_file": "/remote/work/output.out",
            "error_file": None,
            "output_files": (),
            "cache_files": (),
        },
    }

    import asyncio

    runner = asyncio.run(
        context.attach(
            reference,
            _Reader(),
            callback=lambda result: events.append(("finished", None)),
            output_listener=lambda text, kind: events.append((text, kind)),
            request=request,
        )
    )
    runner._done.wait(1)

    assert events == [("remote replay\n", "out"), ("finished", None)]


def test_psij_attach_installs_listener_before_monitor(monkeypatch, tmp_path):
    events = []

    def immediate_monitor(self):
        self._emit_output(b"psij replay\n", "out")
        if self.callback is not None:
            self.callback(self.result)
        self._done.set()

    monkeypatch.setattr(PsijRunner, "_thread_main", immediate_monitor)

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

    reference = RunReference(
        executor="psij:slurm",
        native_id="42",
        directory=str(tmp_path),
        output_file=str(tmp_path / "output.out"),
        error_file=str(tmp_path / "error.err"),
    )
    runner = executor.attach(
        reference,
        _Reader(),
        callback=lambda result: events.append(("finished", None)),
        output_listener=lambda text, kind: events.append((text, kind)),
    )
    runner._done.wait(1)

    assert events == [("psij replay\n", "out"), ("finished", None)]
