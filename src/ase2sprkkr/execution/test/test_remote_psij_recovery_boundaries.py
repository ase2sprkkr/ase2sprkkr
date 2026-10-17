from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.execution.detached import _read_state, _write_state
from ase2sprkkr.execution.remote_psij import _ActiveRun, _RemotePsijAgent


def _state(tmp_path, *, state=RunState.SUBMITTING, **changes):
    run_id = "run"
    run_dir = tmp_path / ".ase2sprkkr" / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    values = {
        "id": run_id,
        "state": state.value,
        "directory": str(tmp_path),
        "backend_reference": None,
        "submitted_at": None,
        "started_at": None,
        "finished_at": None,
        "cancel_requested": False,
    }
    values.update(changes)
    return _write_state(tmp_path, run_id, **values)


def test_executor_construction_failure_during_recovery_is_retryable(tmp_path, monkeypatch):
    state = _state(tmp_path, submitted_at="2026-10-16T12:00:00+00:00")
    agent = _RemotePsijAgent(tmp_path)

    def fail_executor(request):
        raise RuntimeError("PSI/J backend temporarily unavailable")

    monkeypatch.setattr(agent, "_executor_for", fail_executor)

    assert not agent._recover_submitting("run", state, {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.SUBMITTING.value
    assert recovered["submitted_at"] == "2026-10-16T12:00:00+00:00"
    assert recovered["backend_reference"] is None
    assert recovered.get("failure_phase") is None
    assert "temporarily unavailable" in recovered["recovery_error"]
    assert "run" in agent._retry


def test_executor_construction_failure_during_attach_is_retryable(tmp_path, monkeypatch):
    reference = RunReference(
        executor="psij:slurm",
        native_id="job-1",
        directory=str(tmp_path),
    )
    state = _state(
        tmp_path,
        state=RunState.RUNNING,
        backend_reference={
            "executor": reference.executor,
            "native_id": reference.native_id,
            "computer": reference.computer,
            "directory": reference.directory,
            "input_file": reference.input_file,
            "output_file": reference.output_file,
            "error_file": reference.error_file,
        },
        submitted_at="2026-10-16T12:00:00+00:00",
        started_at="2026-10-16T12:01:00+00:00",
    )
    agent = _RemotePsijAgent(tmp_path)

    def fail_executor(request):
        raise RuntimeError("PSI/J backend temporarily unavailable")

    monkeypatch.setattr(agent, "_executor_for", fail_executor)

    assert agent._attach("run", state, {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.RUNNING.value
    assert recovered["finished_at"] is None
    assert "temporarily unavailable" in recovered["attach_error"]
    assert "run" in agent._retry


class _CancelFailingRunner:
    state = RunState.QUEUED

    def cancel(self):
        raise OSError("scheduler cancellation temporarily unavailable")


def test_cancel_backend_error_keeps_active_scheduler_run_retryable(tmp_path):
    _state(
        tmp_path,
        state=RunState.QUEUED,
        cancel_requested=True,
        submitted_at="2026-10-16T12:00:00+00:00",
    )
    agent = _RemotePsijAgent(tmp_path)
    active = _ActiveRun(_CancelFailingRunner(), RunState.QUEUED)
    agent._active["run"] = active

    assert agent._poll_active()

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.QUEUED.value
    assert recovered["finished_at"] is None
    assert "cancellation temporarily unavailable" in recovered["poll_error"]
    assert "run" in agent._active
    assert "run" in agent._retry
    assert not active.cancel_sent


class _PollFailingRunner:
    @property
    def state(self):
        raise OSError("scheduler status temporarily unavailable")


def test_poll_backend_error_keeps_active_scheduler_run_retryable(tmp_path):
    _state(
        tmp_path,
        state=RunState.RUNNING,
        submitted_at="2026-10-16T12:00:00+00:00",
        started_at="2026-10-16T12:01:00+00:00",
    )
    agent = _RemotePsijAgent(tmp_path)
    agent._active["run"] = _ActiveRun(_PollFailingRunner(), RunState.RUNNING)

    assert agent._poll_active()

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.RUNNING.value
    assert recovered["finished_at"] is None
    assert "status temporarily unavailable" in recovered["poll_error"]
    assert "run" in agent._active
    assert "run" in agent._retry
