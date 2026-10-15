from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionManager,
    ExecutionSpec,
    Executor,
    LaunchSpec,
    RunState,
    RunStore,
)
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Runner:
    state = RunState.NEW
    reference = None
    result = None

    def submit(self):
        pass

    def wait(self):
        return None

    def run(self):
        return None

    def cancel(self):
        pass

    def add_output_listener(self, callback):
        pass


class _CapturingExecutor(Executor):
    name = "capture"

    def __init__(self):
        self.run_id = None

    def create_persistent_run(self, reader, spec, *, run_id, callback=None):
        self.run_id = run_id
        return _Runner()

    def create_run(self, reader, spec, callback=None):
        raise AssertionError("persistent manager must use create_persistent_run")


def test_manager_passes_its_identity_to_persistent_executor(tmp_path):
    executor = _CapturingExecutor()
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    spec = SimpleNamespace(directory=str(tmp_path))

    run = manager.create_run(executor, object(), spec, task="DOS", dataset="case")

    assert executor.run_id == run.id


class _Transport:
    host = "cluster"


def test_remote_persistent_runner_uses_manager_run_id(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    reader = SimpleNamespace(create_result=lambda output, input_: object())
    spec = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=str(tmp_path),
        launch=LaunchSpec(local_root=tmp_path, remote_root="/scratch/work"),
    )

    runner = executor.create_persistent_run(
        reader, spec, run_id="ase2sprkkr-run-id"
    )

    assert runner._run_id == "ase2sprkkr-run-id"
