import json
import signal

from ase2sprkkr.execution import FailurePhase, RunState
from ase2sprkkr.execution import remote_service


def _process_handle(tmp_path):
    outcome = tmp_path / "outcome.json"
    handle = {
        "kind": "process",
        "pid": 123,
        "start_time": "456",
        "outcome_file": str(outcome),
        "directory": str(tmp_path),
        "input_file": str(tmp_path / "input"),
        "output_file": str(tmp_path / "output"),
        "error_file": str(tmp_path / "error"),
    }
    return handle, remote_service._encode_handle(handle), outcome


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


def test_status_distinguishes_launch_and_execution_failure(tmp_path, monkeypatch):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)

    _write_outcome(outcome, RunState.FAILED, 127, FailurePhase.LAUNCH)
    launch = remote_service.status_remote(remote_id)
    assert launch["state"] == RunState.FAILED.value
    assert launch["failure_phase"] == FailurePhase.LAUNCH.value

    _write_outcome(outcome, RunState.FAILED, 3, FailurePhase.EXECUTION)
    execution = remote_service.status_remote(remote_id)
    assert execution["state"] == RunState.FAILED.value
    assert execution["return_code"] == 3
    assert execution["failure_phase"] == FailurePhase.EXECUTION.value


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

    cancelled = remote_service.cancel_remote(remote_id)

    assert signals == [(123, signal.SIGTERM)]
    assert cancelled["state"] == RunState.CANCELLED.value
    assert cancelled["failure_phase"] is None

    # A later status query reads the durable backend outcome, rather than
    # regressing to UNKNOWN after the supervisor process is gone.
    assert outcome.is_file()
    later = remote_service.status_remote(remote_id)
    assert later["state"] == RunState.CANCELLED.value
    assert later["failure_phase"] is None
