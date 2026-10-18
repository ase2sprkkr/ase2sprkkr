import io
import json
import signal

import pytest

from ase2sprkkr.execution import FailurePhase, RunState
from ase2sprkkr.execution import psij, remote_process_worker, remote_service


def _process_handle(tmp_path):
    outcome = tmp_path / "remote-outcomes" / ("a" * 32 + ".json")
    handle = {
        "kind": "process",
        "pid": 123,
        "start_time": "456",
        "receipt_id": "a" * 32,
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
    from ase2sprkkr.execution import remote_mailbox
    return remote_mailbox.publish_outcome("a" * 32, state, return_code, failure_phase)


def test_remote_id_is_opaque_backend_handle_without_lifecycle_state(tmp_path):
    handle, remote_id, _ = _process_handle(tmp_path)

    assert remote_id != str(handle)
    decoded = remote_service._decode_handle(remote_id)
    assert decoded["kind"] == "process"
    assert decoded["pid"] == 123
    assert "state" not in decoded
    assert "failure_phase" not in decoded
    assert set(decoded) == {"kind", "pid", "start_time", "receipt_id"}


def test_mailbox_receipt_survives_until_client_ack(tmp_path, monkeypatch):
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)

    assert remote_service.status_remote(remote_id)["state"] == RunState.UNKNOWN.value
    _write_outcome(outcome, RunState.FINISHED, 0)
    assert remote_service.status_remote(remote_id)["state"] == RunState.FINISHED.value
    assert remote_service.status_remote(remote_id)["state"] == RunState.FINISHED.value
    assert remote_service.acknowledge_remote(remote_id) == {"acknowledged": True}
    assert remote_service.acknowledge_remote(remote_id) == {"acknowledged": False}
    assert not outcome.exists()


def test_status_distinguishes_launch_and_execution_failure(tmp_path, monkeypatch):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: False)

    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    _write_outcome(outcome, RunState.FAILED, 127, FailurePhase.LAUNCH)
    launch = remote_service.status_remote(remote_id)
    assert launch["state"] == RunState.FAILED.value
    assert launch["failure_phase"] == FailurePhase.LAUNCH.value

    remote_service.acknowledge_remote(remote_id)
    _write_outcome(outcome, RunState.FAILED, 3, FailurePhase.EXECUTION)
    execution = remote_service.status_remote(remote_id)
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

    def unknown(_remote_id):
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


def test_supervisor_writes_mailbox_receipt(monkeypatch):
    calls = []
    monkeypatch.setattr(
        remote_process_worker,
        "publish_outcome",
        lambda *args: calls.append(args),
    )
    remote_process_worker._write_outcome("a" * 32, RunState.FAILED, 3, FailurePhase.EXECUTION)
    assert calls == [("a" * 32, RunState.FAILED, 3, FailurePhase.EXECUTION)]


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

    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
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


def test_cancel_does_not_publish_terminal_fact_while_process_survives_sigkill(
    tmp_path, monkeypatch
):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setattr(remote_service, "_same_process", lambda _handle: True)
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: True)
    signals = []
    monkeypatch.setattr(
        remote_service.os, "killpg",
        lambda pid, sig: signals.append(sig),
    )

    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    status = remote_service.cancel_remote(remote_id, timeout=0)

    assert status["state"] == RunState.RUNNING.value
    assert signals == [signal.SIGTERM, signal.SIGKILL]
    assert not outcome.exists()


def test_supervisor_requires_parent_launch_release(monkeypatch):
    import os

    read_fd, write_fd = os.pipe()
    os.close(write_fd)
    monkeypatch.setenv(remote_process_worker._START_FD_ENV, str(read_fd))

    with pytest.raises(RuntimeError, match="not released"):
        remote_process_worker._wait_for_parent_release()


def test_uncertain_submit_uses_only_temporary_identity_ticket(tmp_path, monkeypatch):
    from uuid import uuid4
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    run_id = uuid4().hex
    handle = remote_service._encode_handle({
        "kind": "process", "pid": 123,
        "start_time": "456", "receipt_id": "a" * 32,
        "socket_path": "/tmp/fake-socket",
    })
    reply = {"remote_id": handle, "directory": str(tmp_path), "state": "running"}
    remote_service._publish_submission(run_id, reply)
    recovered = remote_service.recover_remote(run_id, {
        "executor": "local", "directory": str(tmp_path),
    })
    assert recovered["remote_id"] == handle
    assert remote_service.acknowledge_submission(run_id) == {"acknowledged": True}
    assert remote_service.recover_remote(run_id, {"executor": "local"}) is None
    assert remote_service.acknowledge_submission(run_id) == {"acknowledged": False}
