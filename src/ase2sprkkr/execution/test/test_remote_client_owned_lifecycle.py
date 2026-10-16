import os
import signal
import sys
import time

import pytest

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.manager import RunRecord, RunStore, SubmissionKey
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor
from ase2sprkkr.execution.remote_service import (
    _decode_handle,
    status_remote,
    submit_remote,
)
from ase2sprkkr.execution.transport import LocalTransport


def _request(tmp_path, code, *, run_id="client-run-id"):
    output = tmp_path / "out"
    return {
        "id": run_id,
        "executor": "local",
        "directory": str(tmp_path),
        "argv": [sys.executable, "-c", code],
        "input_file": None,
        "stdin": None,
        "output_file": str(output),
        "error_file": str(output) + ".err",
        "launch": {},
        "mpi_argv": None,
    }


def _terminal(remote_id):
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        status = status_remote(remote_id)
        if status["state"] != RunState.RUNNING.value:
            return status
        time.sleep(0.02)
    raise AssertionError("remote process did not finish")


def test_remote_process_submission_creates_no_lifecycle_metadata(tmp_path):
    status = submit_remote(
        _request(tmp_path, "import time; time.sleep(0.5)")
    )

    handle = _decode_handle(status["remote_id"])
    try:
        assert handle["kind"] == "process"
        assert handle["pid"] > 0
        assert handle["start_time"]
        assert status_remote(status["remote_id"])["state"] == RunState.RUNNING.value
        assert not (tmp_path / ".ase2sprkkr" / "runs").exists()
    finally:
        try:
            os.killpg(handle["pid"], signal.SIGKILL)
        except ProcessLookupError:
            pass


def test_remote_process_preserves_success_return_code(tmp_path):
    submitted = submit_remote(_request(tmp_path, "print('ok')"))

    status = _terminal(submitted["remote_id"])

    assert status["state"] == RunState.FINISHED.value
    assert status["return_code"] == 0
    assert not (tmp_path / ".ase2sprkkr" / "runs").exists()


def test_remote_process_preserves_failure_return_code(tmp_path):
    submitted = submit_remote(_request(tmp_path, "raise SystemExit(7)"))

    status = _terminal(submitted["remote_id"])

    assert status["state"] == RunState.EXECUTION_FAILED.value
    assert status["return_code"] == 7


def test_remote_process_reports_child_launch_failure(tmp_path):
    request = _request(tmp_path, "unused")
    request["argv"] = [str(tmp_path / "does-not-exist")]

    submitted = submit_remote(request)
    status = _terminal(submitted["remote_id"])

    assert status["state"] == RunState.EXECUTION_FAILED.value
    assert status["return_code"] == 127


def test_run_reference_separates_remote_id_from_client_context():
    reference = RunReference(
        executor="remote:local",
        native_id="opaque-target-id",
        context={
            "launch": {"remote_root": "/scratch/user"},
            "client_sync": {"output_file": "/scratch/user/case/out"},
        },
    )

    assert reference.native_id == "opaque-target-id"
    assert reference.context == {
        "launch": {"remote_root": "/scratch/user"},
        "client_sync": {"output_file": "/scratch/user/case/out"},
    }


def test_client_reconnect_context_roundtrips_through_run_store(tmp_path):
    endpoint = ExecutorEndpoint(executor="ssh", host="cluster")
    submission = SubmissionKey(
        endpoint=endpoint,
        directory="/work/case",
        task="SCF",
        context={"launch": {"remote_root": "/scratch/user"}},
    )
    reference = RunReference(
        executor="remote:local",
        native_id="opaque-target-id",
        computer="cluster",
        directory="/work/case",
        context={"client_sync": {"output_file": "/scratch/user/case/out"}},
    )
    record = RunRecord(
        id="run-1",
        state=RunState.RUNNING,
        created_at="2026-10-17T00:00:00+00:00",
        updated_at="2026-10-17T00:00:00+00:00",
        task="SCF",
        directory=str(tmp_path),
        endpoint=endpoint,
        submission=submission,
        reference=reference,
    )
    store = RunStore(tmp_path / "runs")

    store.save(record)
    loaded = store.get("run-1")

    assert loaded.reference.native_id == "opaque-target-id"
    assert loaded.reference.context == reference.context
    assert loaded.submission.context == submission.context


def test_remote_endpoint_status_sends_only_opaque_remote_id():
    class Transport:
        host = "cluster"

    executor = RemoteEndpointExecutor(Transport())
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {"state": RunState.RUNNING.value, "return_code": None}

    executor._json_request = request
    reference = RunReference(
        executor="remote:local",
        native_id="pid:start",
        computer="cluster",
        directory="/work/case",
        context={"launch": {}, "client_sync": {}},
    )
    status = executor.status(reference)[0]

    assert status.state == RunState.RUNNING
    assert calls == [("status", ("pid:start",), None)]


def test_direct_remote_submission_recovery_is_explicitly_unsupported(tmp_path):
    executor = RemoteEndpointExecutor(LocalTransport())

    class Submission:
        directory = str(tmp_path)
        context = {"launch": {}}

    async def request(action, *args, payload=None):
        assert action == "recover"
        return {
            "unsupported": True,
            "error": "bare process recovery is unavailable",
        }

    executor._json_request = request
    with pytest.raises(NotImplementedError, match="bare process recovery"):
        executor.recover_submission("run", Submission())
