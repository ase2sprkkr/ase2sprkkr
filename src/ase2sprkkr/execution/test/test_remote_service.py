import io
import json

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
        "socket_path": "/tmp/test-supervisor.sock",
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
    assert set(decoded) == {"kind", "pid", "start_time", "receipt_id", "socket_path"}


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
        remote_service, "_psij_executor",
        lambda opts: pytest.fail("local Slurm status must not initialize PSI/J"),
    )
    monkeypatch.setattr(
        psij,
        "_slurm_details",
        lambda native_id: {"JobState": "RUNNING", "ExitCode": "0:0"},
    )

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.RUNNING.value
    assert status["return_code"] is None
    assert status["failure_phase"] is None


def test_slurm_status_falls_back_to_exact_accounting_record(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(psij, "_slurm_details", lambda _native_id: {})
    monkeypatch.setattr(psij, "_scheduler_output", lambda argv: "42|COMPLETED|0:0|\n")

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.FINISHED.value
    assert status["return_code"] == 0
    assert status["failure_phase"] is None


def test_slurm_failed_accounting_record_is_execution_failure(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(psij, "_slurm_details", lambda _native_id: {})
    monkeypatch.setattr(psij, "_scheduler_output", lambda argv: "42|FAILED|7:0|\n")

    status = remote_service.status_remote(remote_id)

    assert status["state"] == RunState.FAILED.value
    assert status["return_code"] == 7
    assert status["failure_phase"] == FailurePhase.EXECUTION.value



def test_slurm_completed_with_nonzero_exit_is_failed(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    monkeypatch.setattr(psij, "_slurm_details", lambda _id: {})
    monkeypatch.setattr(
        psij, "_scheduler_output", lambda argv: "42|COMPLETED|7:0|\n",
    )
    status = remote_service.status_remote(remote_id)
    assert status["state"] == RunState.FAILED.value
    assert status["return_code"] == 7
    assert status["failure_phase"] == FailurePhase.EXECUTION.value


def test_psij_cancellation_does_not_invent_terminal_state(monkeypatch, tmp_path):
    _, remote_id = _psij_handle(tmp_path)
    calls = []

    class Job:
        def cancel(self):
            calls.append("cancel")

    class Backend:
        _psij = type("Psij", (), {"Job": Job})
        _job_executor = type("Executor", (), {
            "attach": lambda self, job, native_id: calls.append(native_id)
        })()

    monkeypatch.setattr(remote_service, "_psij_executor", lambda _: Backend())
    monkeypatch.setattr(remote_service, "status_remote", lambda _: {
        "state": RunState.RUNNING.value, "return_code": None,
        "failure_phase": None,
    })
    status = remote_service.cancel_remote(remote_id)
    assert status["state"] == RunState.RUNNING.value
    assert calls == ["42", "cancel"]


def test_follow_does_not_interpret_unknown_as_terminal(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    states = iter((
        RunState.UNKNOWN.value, RunState.UNKNOWN.value,
        RunState.FINISHED.value,
    ))
    seen = []
    def status(_remote_id):
        seen.append(_remote_id)
        return {"state": next(states)}
    monkeypatch.setattr(remote_service, "status_remote", status)
    outcome = remote_service.follow_remote(
        remote_id, io.BytesIO(), io.BytesIO(), poll_interval=0.001,
    )
    assert outcome["state"] == RunState.FINISHED.value
    assert seen == [remote_id] * 3


def test_follow_transient_unknown_resolves_before_terminal_state(tmp_path, monkeypatch):
    _, remote_id = _psij_handle(tmp_path)
    states = iter((
        RunState.UNKNOWN.value, RunState.RUNNING.value,
        RunState.FINISHED.value,
    ))
    seen = []
    def status(_id):
        seen.append(_id)
        return {"state": next(states)}
    monkeypatch.setattr(remote_service, "status_remote", status)
    monkeypatch.setattr(remote_service.time, "sleep", lambda _: None)
    outcome = remote_service.follow_remote(
        remote_id, io.BytesIO(), io.BytesIO(), poll_interval=0,
    )
    assert outcome["state"] == RunState.FINISHED.value
    assert seen == [remote_id] * 3


def test_supervisor_writes_mailbox_receipt(monkeypatch):
    calls = []
    monkeypatch.setattr(
        remote_process_worker,
        "publish_outcome",
        lambda *args: calls.append(args),
    )
    remote_process_worker._write_outcome("a" * 32, RunState.FAILED, 3, FailurePhase.EXECUTION)
    assert calls == [("a" * 32, RunState.FAILED, 3, FailurePhase.EXECUTION)]


def test_cancel_requests_supervisor_without_inventing_terminal_state(
    tmp_path, monkeypatch
):
    _, remote_id, outcome = _process_handle(tmp_path)
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(remote_service, "_process_running", lambda _handle: True)
    sent = []

    class Connection:
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def settimeout(self, timeout):
            pass
        def connect(self, path):
            assert path == "/tmp/test-supervisor.sock"
        def sendall(self, data):
            sent.append(data)
        def recv(self, count):
            return b'{"accepted":true}\n'

    monkeypatch.setattr(remote_service.socket, "socket", lambda *args: Connection())
    status = remote_service.cancel_remote(remote_id)
    assert sent == [b'{"command":"cancel"}\n']
    assert status["state"] == RunState.RUNNING.value
    assert not outcome.exists()


def test_cancel_after_natural_completion_keeps_backend_outcome(tmp_path, monkeypatch):
    _, remote_id, _ = _process_handle(tmp_path)
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    _write_outcome(None, RunState.FINISHED, 0)

    def forbidden(*args):
        raise AssertionError("Terminal process must not receive cancellation")
    monkeypatch.setattr(remote_service.socket, "socket", forbidden)

    status = remote_service.cancel_remote(remote_id)
    assert status["state"] == RunState.FINISHED.value
    assert status["return_code"] == 0


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


def test_direct_follow_uses_shared_socket_protocol(tmp_path, monkeypatch):
    """The target command and the supervisor must agree on frames and offsets."""
    import asyncio
    import os
    import sys

    async def exercise():
        monkeypatch.setattr(
            remote_process_worker, "publish_outcome", lambda *args, **kwargs: None,
        )
        rd, wr = os.pipe()
        spec = {
            "argv": [sys.executable, "-c",
                     "import sys,time; print('stdout', flush=True); "
                     "print('stderr', file=sys.stderr, flush=True); time.sleep(.1)"],
            "directory": str(tmp_path),
            "output_file": str(tmp_path / "stdout"),
            "error_file": str(tmp_path / "stderr"),
            "receipt_id": "a" * 32,
            "socket_fd": wr,
            "detached": True,
        }
        worker = asyncio.create_task(remote_process_worker._run(spec))
        with os.fdopen(rd, "rb") as source:
            socket_path = (await asyncio.to_thread(source.readline)).decode().strip()
        remote_id = remote_service._encode_handle({
            "kind": "process", "socket_path": socket_path,
            "pid": os.getpid(), "start_time": "0", "receipt_id": "a" * 32,
        })
        out, err = io.BytesIO(), io.BytesIO()
        terminal = await asyncio.to_thread(
            remote_service.follow_remote, remote_id, out, err,
            stdout_path=spec["output_file"], stderr_path=spec["error_file"],
        )
        assert terminal["state"] == "finished"
        assert out.getvalue() == b"stdout\\n"
        assert err.getvalue() == b"stderr\\n"
        assert await worker == 0

    asyncio.run(exercise())


def test_direct_follow_without_stream_files_is_live_only(tmp_path, monkeypatch):
    import asyncio
    import os
    import sys

    async def exercise():
        monkeypatch.setattr(
            remote_process_worker, "publish_outcome", lambda *args, **kwargs: None,
        )
        rd, wr = os.pipe()
        spec = {
            "argv": [sys.executable, "-c",
                     "import time; print('live', flush=True); time.sleep(.2)"],
            "directory": str(tmp_path),
            "output_file": None, "error_file": None,
            "receipt_id": "b" * 32,
            "socket_fd": wr, "detached": True,
        }
        worker = asyncio.create_task(remote_process_worker._run(spec))
        with os.fdopen(rd, "rb") as source:
            socket_path = (await asyncio.to_thread(source.readline)).decode().strip()
        remote_id = remote_service._encode_handle({
            "kind": "process", "socket_path": socket_path,
            "pid": os.getpid(), "start_time": "0", "receipt_id": "b" * 32,
        })
        out, err = io.BytesIO(), io.BytesIO()
        terminal = await asyncio.to_thread(
            remote_service.follow_remote, remote_id, out, err,
            stdout_offset=0, stderr_offset=0,
        )
        assert terminal["state"] == "finished"
        assert out.getvalue() == b"live\\n"
        assert err.getvalue() == b""
        assert await worker == 0
        assert list(tmp_path.iterdir()) == []

    asyncio.run(exercise())
