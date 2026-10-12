import json
from types import SimpleNamespace

from ..base import RunReference, RunState
from ..remote import _request_path, _run_directory, _write_state, status_remote
from ..remote_psij import (
    _RemotePsijAgent,
    cancel_remote_psij,
    ensure_remote_psij_agent,
)


def _request(workspace):
    return {
        "executor": "psij",
        "executor_options": {"backend": "slurm", "url": None},
        "directory": str(workspace),
        "argv": ["true"],
        "input_file": None,
        "stdin": None,
        "output_file": str(workspace / "out"),
        "error_file": str(workspace / "err"),
    }


def _record(workspace, run_id, *, state="queued", backend_reference=None):
    run_dir = _run_directory(workspace, run_id)
    run_dir.mkdir(parents=True)
    request = dict(_request(workspace), id=run_id)
    _request_path(workspace, run_id).write_text(json.dumps(request))
    _write_state(
        workspace,
        run_id,
        id=run_id,
        state=state,
        directory=str(workspace),
        pid=None,
        backend_reference=backend_reference,
        finished_at=None,
    )
    return request


class _FakeRunner:
    def __init__(self, workspace, native_id, *, state=RunState.QUEUED):
        self._state = state
        self.result = SimpleNamespace(return_code=0)
        self.reference = RunReference(
            executor="psij:slurm",
            native_id=native_id,
            directory=str(workspace),
            output_file=str(workspace / (native_id + ".out")),
            error_file=str(workspace / (native_id + ".err")),
        )
        self.submit_count = 0
        self.cancel_count = 0

    @property
    def state(self):
        return self._state

    def submit(self):
        self.submit_count += 1

    def cancel(self):
        self.cancel_count += 1
        self._state = RunState.CANCELLED

    def wait(self):
        return self.result


class _FakeExecutor:
    def __init__(self, workspace):
        self.workspace = workspace
        self.created = []
        self.attached = []

    def create_run(self, reader, spec):
        runner = _FakeRunner(self.workspace, "job-{}".format(len(self.created) + 1))
        self.created.append(runner)
        return runner

    def attach(self, reference, reader):
        runner = _FakeRunner(
            self.workspace,
            reference.native_id,
            state=RunState.RUNNING,
        )
        self.attached.append((reference, runner))
        return runner


def test_remote_psij_agent_reuses_executor_for_multiple_runs(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    _record(workspace, "run-1")
    _record(workspace, "run-2")

    executor = _FakeExecutor(workspace)
    created_executors = []

    def factory(request):
        created_executors.append(request["executor_options"])
        return executor

    monkeypatch.setattr("ase2sprkkr.execution.remote_psij._executor", factory)

    agent = _RemotePsijAgent(workspace)
    assert agent._discover()

    assert created_executors == [{"backend": "slurm", "url": None}]
    assert len(executor.created) == 2
    assert status_remote(workspace, "run-1")["backend_reference"]["native_id"] == "job-1"
    assert status_remote(workspace, "run-2")["backend_reference"]["native_id"] == "job-2"


def test_remote_psij_agent_recovers_existing_scheduler_job(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    reference = RunReference(
        executor="psij:slurm",
        native_id="job-existing",
        directory=str(workspace),
        output_file=str(workspace / "out"),
        error_file=str(workspace / "err"),
    )
    _record(
        workspace,
        "run-1",
        state="running",
        backend_reference={
            "executor": reference.executor,
            "native_id": reference.native_id,
            "computer": reference.computer,
            "directory": reference.directory,
            "input_file": reference.input_file,
            "output_file": reference.output_file,
            "error_file": reference.error_file,
        },
    )

    executor = _FakeExecutor(workspace)
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._executor", lambda request: executor
    )

    agent = _RemotePsijAgent(workspace)
    assert agent._discover()

    assert executor.created == []
    assert len(executor.attached) == 1
    assert executor.attached[0][0].native_id == "job-existing"


def test_remote_psij_cancel_is_serviced_by_shared_agent(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    request = _record(workspace, "run-1")

    executor = _FakeExecutor(workspace)
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._executor", lambda request: executor
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij.ensure_remote_psij_agent",
        lambda directory, run_ids: True,
    )

    agent = _RemotePsijAgent(workspace)
    assert agent._discover()
    runner = executor.created[0]

    state = cancel_remote_psij(workspace, "run-1", request, timeout=0)
    assert state["cancel_requested"] is True
    assert runner.cancel_count == 0

    agent._poll_active()

    assert runner.cancel_count == 1
    assert status_remote(workspace, "run-1")["state"] == "cancelled"


def test_remote_psij_status_recovery_starts_agent_for_active_run(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    _record(workspace, "run-1", state="running")
    spawned = []

    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._agent_is_running", lambda directory: False
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._spawn_agent",
        lambda directory, run_id: spawned.append((directory, run_id)),
    )

    assert ensure_remote_psij_agent(workspace, ("run-1",))
    assert spawned == [(workspace.resolve(), "run-1")]
