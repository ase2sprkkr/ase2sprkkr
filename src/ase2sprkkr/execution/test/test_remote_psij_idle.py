import json

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.execution.remote import _request_path, _run_directory, _write_state
from ase2sprkkr.execution.remote_psij import _RemotePsijAgent


class _RunningRunner:
    def __init__(self, reference):
        self.reference = reference
        self.state = RunState.RUNNING

    def cancel(self):
        self.state = RunState.CANCELLED

    def wait(self):
        raise AssertionError("running scheduler job must not be waited during idle exit")


class _Executor:
    def __init__(self):
        self.attached = []

    def attach(self, reference, reader):
        self.attached.append(reference)
        return _RunningRunner(reference)


def _record(workspace):
    run_id = "run-1"
    run_dir = _run_directory(workspace, run_id)
    run_dir.mkdir(parents=True)
    request = {
        "id": run_id,
        "executor": "psij",
        "executor_options": {"backend": "slurm", "url": None},
        "directory": str(workspace),
        "argv": ["true"],
        "input_file": None,
        "stdin": None,
        "output_file": str(workspace / "out"),
        "error_file": str(workspace / "err"),
    }
    _request_path(workspace, run_id).write_text(json.dumps(request))
    reference = RunReference(
        executor="psij:slurm",
        native_id="job-1",
        directory=str(workspace),
        output_file=str(workspace / "out"),
        error_file=str(workspace / "err"),
    )
    _write_state(
        workspace,
        run_id,
        id=run_id,
        state=RunState.RUNNING.value,
        directory=str(workspace),
        pid=None,
        backend_reference={
            "executor": reference.executor,
            "native_id": reference.native_id,
            "computer": reference.computer,
            "directory": reference.directory,
            "input_file": reference.input_file,
            "output_file": reference.output_file,
            "error_file": reference.error_file,
            "target": reference.target,
        },
        created_at="2026-10-12T00:00:00+00:00",
        finished_at=None,
    )
    return run_id


def test_agent_exits_when_scheduler_job_is_still_running(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    run_id = _record(workspace)
    executor = _Executor()
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._executor",
        lambda request: executor,
    )
    monkeypatch.setattr("ase2sprkkr.execution.remote_psij._AGENT_IDLE_TIMEOUT", 0.0)
    monkeypatch.setattr("ase2sprkkr.execution.remote_psij._AGENT_POLL_INTERVAL", 0.0)

    agent = _RemotePsijAgent(workspace)
    agent.run()

    assert [reference.native_id for reference in executor.attached] == ["job-1"]
    assert run_id in agent._active
