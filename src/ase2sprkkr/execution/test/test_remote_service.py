import io
import json
import signal

from ase2sprkkr.execution import FailurePhase, RunState
from ase2sprkkr.execution import psij, remote_process_worker, remote_service


def _process_handle(tmp_path):
    outcome = tmp_path / "outcome.json"
    handle = {
        "kind": "process",
        "pid": 123,
        "start_time": "456",
    }
    return handle, remote_service._encode_handle(handle), outcome


def _psij_handle(tmp_path):
    handle = {
        "kind": "psij",
        "backend": "slurm",
        "url": None,
        "native_id": "42",
    }
    return handle, remote_service._encode_handle(handle)


def _write_outcome(path, state, return_code, failure_phase=None):
    path.write_text(
        json.dumps(
            {
                "state": state.value,
                "return_code": return_code,
                "failure_phase": (
                    None if failure_phase is None else failure_phase.value
                ),
            }
        )
    )


def test_remote_id_is_opaque_backend_handle_without_lifecycle_state(tmp_path):
    handle, remote_id, _ = _process_handle(tmp_path)

    assert remote_id != str(handle)
    decoded = remote_service._decode_handle(remote_id)
    assert decoded["kind"] == "process"
    assert decoded["pid"] == 123
    assert "state" not in decoded
    assert "failure_phase" not in decoded
    assert set(decoded) == {"kind", "pid", "start_time"}


def test_direct_process_without_client_receipt_never_infers_success(tmp_path, monkeypatch):
    _, remote_id, outcome = _process_handle(tmp_path)
    _write_outcome(outcome, RunState.FINISHED, 0)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)

    assert remote_service.status_remote(remote_id)["state"] == RunState.UNKNOWN.value
    assert remote_service.status_remote(remote_id, outcome_file=outcome) == {
        "state": RunState.FINISHED.value,
        "return_code": 0,
        "failure_phase": None,
    }


def test_status_distinguishes_launch_and_execution_failure(tmp_path, monkeypatch):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)

    _write_outcome(outcome, RunState.FAILED, 127, FailurePhase.LAUNCH)
    launch = remote_service.status_remote(remote_id, outcome_file=outcome)
    assert launch["state"] == RunState.FAILED.value
    assert launch["failure_phase"] == FailurePhase.LAUNCH.value

    _write_outcome(outcome, RunState.FAILED, 3, FailurePhase.EXECUTION)
    execution = remote_service.status_remote(remote_id, outcome_file=outcome)
    assert execution["state"] == RunState.FAILED.value
    assert execution["return_code"] == 3
    assert execution["failure_phase"] == FailurePhase.EXECUTION.value


def test_slurm_status_is_native_and_does_not_use_psij_job_cache(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(
        psij,
        "_slurm_details",
        lambda native_id: {"JobState": "RUNNING", "ExitCode": "0:0"},
    )

    def forbidden_executor(_options):
        raise AssertionError("stateless Slurm status must not create a PSI/J Job")

    monkeypatch.setattr(remote_service, "_psij_executor", forbidden_executor)

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.RUNNING.value
    assert status["return_code"] is None
    assert status["failure_phase"] is None


def test_slurm_status_falls_back_to_exact_accounting_record(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(psij, "_slurm_details", lambda _native_id: {})
    monkeypatch.setattr(
        remote_service,
        "_slurm_accounting_status",
        lambda native_id: (RunState.FINISHED, 0),
    )

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.FINISHED.value
    assert status["return_code"] == 0
    assert status["failure_phase"] is None


def test_slurm_failed_accounting_record_is_execution_failure(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(psij, "_slurm_details", lambda _native_id: {})
    monkeypatch.setattr(
        remote_service,
        "_slurm_accounting_status",
        lambda native_id: (RunState.FAILED, 7),
    )

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.FAILED.value
    assert status["return_code"] == 7
    assert status["failure_phase"] == FailurePhase.EXECUTION.value


def test_follow_returns_unknown_instead_of_waiting_forever(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    calls = []

    def unknown(_remote_id, *, outcome_file=None):
        calls.append(_remote_id)
        return {"state": RunState.UNKNOWN.value}

    monkeypatch.setattr(remote_service, "status_remote", unknown)

    status = remote_service.follow_remote(
        remote_id,
        io.BytesIO(),
        io.BytesIO(),
        poll_interval=0,
    )

    assert status["state"] == RunState.UNKNOWN.value
    assert calls == [remote_id]


def test_direct_process_outcome_uses_crash_durable_write(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        remote_process_worker,
        "durable_write_json",
        lambda path, payload, **options: calls.append((path, payload, options)),
    )
    outcome = tmp_path / "outcome.json"

    remote_process_worker._write_outcome(
        outcome,
        RunState.FAILED,
        3,
        FailurePhase.EXECUTION,
    )

    assert calls == [
        (
            outcome,
            {
                "state": RunState.FAILED.value,
                "return_code": 3,
                "failure_phase": FailurePhase.EXECUTION.value,
            },
            {"separators": (",", ":")},
        )
    ]


def test_cancel_publishes_durable_cancelled_backend_outcome(
    tmp_path, monkeypatch
):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_same_process", lambda _handle: True)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)
    signals = []
    monkeypatch.setattr(
        remote_service.os,
        "killpg",
        lambda pid, sig: signals.append((pid, sig)),
    )

    cancelled = remote_service.cancel_remote(remote_id, outcome_file=outcome)

    assert signals == [(123, signal.SIGTERM)]
    assert cancelled["state"] == RunState.CANCELLED.value
    assert cancelled["failure_phase"] is None

    # A later status query reads the durable backend outcome, rather than
    # regressing to UNKNOWN after the supervisor process is gone.
    assert outcome.is_file()
    later = remote_service.status_remote(remote_id, outcome_file=outcome)
    assert later["state"] == RunState.CANCELLED.value
    assert later["failure_phase"] is None
