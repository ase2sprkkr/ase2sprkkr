import sys
import time
from pathlib import Path

import pytest

from ase2sprkkr.execution import ExecutionSpec, LaunchSpec, RunState
from ase2sprkkr.execution.remote_endpoint_client import (
    RemoteEndpointExecutor,
    _pack_reference,
    _unpack_reference,
)
from ase2sprkkr.execution.remote_service import (
    _decode_handle,
    status_remote,
    submit_remote,
)
from ase2sprkkr.execution.transport import LocalTransport


def test_remote_process_submission_creates_no_lifecycle_metadata(tmp_path):
    output = tmp_path / "out"
    status = submit_remote(
        {
            "id": "client-run-id",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [sys.executable, "-c", "import time; time.sleep(0.5)"],
            "input_file": None,
            "stdin": None,
            "output_file": str(output),
            "error_file": str(output) + ".err",
            "launch": {},
            "mpi_argv": None,
        }
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
            import os
            import signal

            os.killpg(handle["pid"], signal.SIGKILL)
        except ProcessLookupError:
            pass


def test_remote_reference_keeps_reconnect_context_on_client():
    native_id = _pack_reference(
        "opaque-target-id",
        launch={"remote_root": "/scratch/user"},
        client_sync={"output_file": "/scratch/user/case/out"},
    )

    assert _unpack_reference(native_id) == {
        "remote_id": "opaque-target-id",
        "launch": {"remote_root": "/scratch/user"},
        "client_sync": {"output_file": "/scratch/user/case/out"},
    }


def test_remote_endpoint_status_sends_only_opaque_remote_id():
    class Transport:
        host = "cluster"

    executor = RemoteEndpointExecutor(Transport())
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {"state": RunState.RUNNING.value, "return_code": None}

    executor._json_request = request
    native_id = _pack_reference("pid:start", launch={}, client_sync={})

    from ase2sprkkr.execution import RunReference

    reference = RunReference(
        executor="remote:local",
        native_id=native_id,
        computer="cluster",
        directory="/work/case",
    )
    status = executor.status(reference)[0]

    assert status.state == RunState.RUNNING
    assert calls == [("status", ("pid:start",), None)]


def test_direct_remote_submission_recovery_is_explicitly_unsupported(tmp_path):
    executor = RemoteEndpointExecutor(LocalTransport())

    class Submission:
        directory = str(tmp_path)

    async def request(action, *args, payload=None):
        assert action == "recover"
        return {
            "unsupported": True,
            "error": "bare process recovery is unavailable",
        }

    executor._json_request = request
    with pytest.raises(NotImplementedError, match="bare process recovery"):
        executor.recover_submission("run", Submission())
