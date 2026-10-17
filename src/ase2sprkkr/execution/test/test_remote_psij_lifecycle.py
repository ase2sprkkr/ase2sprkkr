from types import SimpleNamespace

import ase2sprkkr.execution.detached as detached
from ase2sprkkr.execution import FailurePhase, RunReference, RunState
from ase2sprkkr.execution.detached import _read_state, _write_state
from ase2sprkkr.execution.remote_psij import _RemotePsijAgent


def _initial_state(tmp_path, run_id="run"):
    run_dir = tmp_path / ".ase2sprkkr" / "runs" / run_id
    run_dir.mkdir(parents=True)
    return _write_state(
        tmp_path,
        run_id,
        id=run_id,
        state=RunState.SUBMITTING.value,
        directory=str(tmp_path),
        backend_reference=None,
        submitted_at=None,
        started_at=None,
        finished_at=None,
        cancel_requested=False,
    )


class _Runner:
    def __init__(self, state=RunState.QUEUED):
        self._state = state
        self.reference = RunReference(
            executor="psij:slurm",
            native_id="123",
            directory="/work",
        )

    @property
    def state(self):
        return self._state

    def submit(self):
        return None


class _Executor:
    def __init__(self, runner):
        self.runner = runner

    def create_persistent_run(self, reader, spec, *, run_id):
        return self.runner


def test_psij_becomes_queued_only_after_scheduler_acceptance(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    runner = _Runner(RunState.QUEUED)
    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: _Executor(runner))

    assert agent._submit("run", {"argv": ["true"]})

    state = _read_state(tmp_path, "run")
    assert state["state"] == RunState.QUEUED.value
    assert state["submitted_at"] is not None
    assert state["backend_reference"]["native_id"] == "123"
    assert state["started_at"] is None
    assert state["failure_phase"] is None

    runner._state = RunState.RUNNING
    assert agent._poll_active()

    state = _read_state(tmp_path, "run")
    assert state["state"] == RunState.RUNNING.value
    assert state["started_at"] is not None


def test_psij_backend_reference_write_failure_stays_recoverable(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    runner = _Runner(RunState.QUEUED)
    agent = _RemotePsijAgent(tmp_path)
    executor = _Executor(runner)
    monkeypatch.setattr(agent, "_executor_for", lambda request: executor)

    real_write = detached._write_state_unlocked
    failed = False

    def fail_reference_publication(workspace, run_id, **changes):
        nonlocal failed
        if not failed and changes.get("backend_reference") is not None:
            failed = True
            raise OSError("backend reference write failed")
        return real_write(workspace, run_id, **changes)

    monkeypatch.setattr(detached, "_write_state_unlocked", fail_reference_publication)

    assert agent._submit("run", {"argv": ["true"]})

    state = _read_state(tmp_path, "run")
    assert state["state"] == RunState.SUBMITTING.value
    assert state["submitted_at"] is not None
    assert state["backend_reference"] is None
    assert state.get("failure_phase") is None
    assert "run" in agent._retry

    class RecoveringExecutor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            assert submitted_at is not None
            return SimpleNamespace(
                reference=runner.reference,
                state=RunState.QUEUED,
                return_code=None,
                message=None,
                failure_phase=None,
            )

    monkeypatch.setattr(
        agent, "_executor_for", lambda request: RecoveringExecutor()
    )
    assert agent._recover_submitting("run", state, {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.QUEUED.value
    assert recovered["backend_reference"]["native_id"] == "123"
    assert recovered["failure_phase"] is None


def test_psij_missing_reference_after_submit_stays_recoverable(tmp_path, monkeypatch):
    _initial_state(tmp_path)

    class MissingReferenceRunner(_Runner):
        def __init__(self):
            super().__init__(RunState.QUEUED)
            self.reference = None

    runner = MissingReferenceRunner()
    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: _Executor(runner))

    assert agent._submit("run", {"argv": ["true"]})

    state = _read_state(tmp_path, "run")
    assert state["state"] == RunState.SUBMITTING.value
    assert state["submitted_at"] is not None
    assert state["backend_reference"] is None
    assert state.get("submission_error") is None
    assert "persistent job reference" in state["recovery_error"]
    assert state.get("failure_phase") is None
    assert "run" in agent._retry

    class AbsentExecutor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            assert submitted_at is not None
            return None

    monkeypatch.setattr(agent, "_executor_for", lambda request: AbsentExecutor())
    assert agent._recover_submitting("run", state, {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.SUBMITTING.value
    assert recovered["submitted_at"] is None
    assert recovered["backend_reference"] is None
    assert recovered["finished_at"] is None
    assert recovered["recovery_error"] is None
    assert recovered["error"] is None


def test_psij_incomplete_recovery_stays_submitting(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    _write_state(tmp_path, "run", submitted_at="2026-10-16T12:00:00+00:00")

    class RecoveringExecutor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            return SimpleNamespace(
                reference=None,
                state=RunState.QUEUED,
                return_code=None,
                message=None,
                failure_phase=None,
            )

    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: RecoveringExecutor())

    assert not agent._recover_submitting("run", _read_state(tmp_path, "run"), {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.SUBMITTING.value
    assert recovered["backend_reference"] is None
    assert recovered["finished_at"] is None
    assert "lacks a job reference" in recovered["recovery_error"]
    assert "run" in agent._retry


def test_psij_recovered_reference_write_failure_stays_recoverable(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    _write_state(tmp_path, "run", submitted_at="2026-10-16T12:00:00+00:00")
    reference = _Runner().reference

    class RecoveringExecutor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            assert submitted_at == "2026-10-16T12:00:00+00:00"
            return SimpleNamespace(
                reference=reference,
                state=RunState.QUEUED,
                return_code=None,
                message=None,
                failure_phase=None,
            )

        def attach(self, reference, reader):
            raise OSError("attach deferred")

    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: RecoveringExecutor())

    real_write = detached._write_state_unlocked
    failed = False

    def fail_recovered_reference_publication(workspace, run_id, **changes):
        nonlocal failed
        if not failed and changes.get("backend_reference") is not None:
            failed = True
            raise OSError("recovered reference write failed")
        return real_write(workspace, run_id, **changes)

    monkeypatch.setattr(
        detached,
        "_write_state_unlocked",
        fail_recovered_reference_publication,
    )

    assert not agent._recover_submitting("run", _read_state(tmp_path, "run"), {})

    retryable = _read_state(tmp_path, "run")
    assert retryable["state"] == RunState.SUBMITTING.value
    assert retryable["backend_reference"] is None
    assert retryable.get("failure_phase") is None
    assert "recovered reference write failed" in retryable["recovery_error"]
    assert "run" in agent._retry

    assert agent._recover_submitting("run", retryable, {})
    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.QUEUED.value
    assert recovered["backend_reference"]["native_id"] == "123"
    assert recovered["failure_phase"] is None


def test_psij_proven_absent_handoff_resets_for_resubmission(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    _write_state(
        tmp_path,
        "run",
        submitted_at="2026-10-16T12:00:00+00:00",
        recovery_error="backend reference write failed",
        error="backend reference write failed",
    )

    class Executor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            assert submitted_at == "2026-10-16T12:00:00+00:00"
            return None

    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: Executor())

    assert agent._recover_submitting("run", _read_state(tmp_path, "run"), {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.SUBMITTING.value
    assert recovered["submitted_at"] is None
    assert recovered["backend_reference"] is None
    assert recovered["recovery_error"] is None
    assert recovered["error"] is None
    assert recovered["finished_at"] is None


def test_psij_retry_diagnostic_write_failure_is_not_terminal(tmp_path, monkeypatch):
    _initial_state(tmp_path)

    class AmbiguousRunner(_Runner):
        def submit(self):
            raise OSError("scheduler submit connection dropped")

    runner = AmbiguousRunner()
    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: _Executor(runner))

    real_write = detached._write_state_unlocked

    def fail_diagnostic(workspace, run_id, **changes):
        if changes.get("submission_error") is not None:
            raise OSError("diagnostic write failed")
        return real_write(workspace, run_id, **changes)

    monkeypatch.setattr(detached, "_write_state_unlocked", fail_diagnostic)

    assert agent._submit("run", {"argv": ["true"]})

    state = _read_state(tmp_path, "run")
    assert state["state"] == RunState.SUBMITTING.value
    assert state["submitted_at"] is not None
    assert state["backend_reference"] is None
    assert state.get("failure_phase") is None
    assert "run" in agent._retry


def test_psij_proven_absent_submission_is_submission_failure(tmp_path, monkeypatch):
    _initial_state(tmp_path)
    _write_state(
        tmp_path,
        "run",
        submitted_at="2026-10-16T12:00:00+00:00",
        submission_error="scheduler submit failed",
    )

    class Executor:
        def recover_submission(self, run_id, submission, *, submitted_at=None):
            assert submitted_at == "2026-10-16T12:00:00+00:00"
            return None

    agent = _RemotePsijAgent(tmp_path)
    monkeypatch.setattr(agent, "_executor_for", lambda request: Executor())

    assert agent._recover_submitting("run", _read_state(tmp_path, "run"), {})

    recovered = _read_state(tmp_path, "run")
    assert recovered["state"] == RunState.FAILED.value
    assert recovered["failure_phase"] == FailurePhase.SUBMISSION.value
    assert recovered["started_at"] is None
