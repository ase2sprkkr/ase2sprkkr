import pytest

from ase2sprkkr.execution import (
    Execution, ExecutorEndpoint, FileManifest, RunRecord, FailurePhase, RunState,
    WorkingDirectory,
)
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.execution.remote_endpoint import (
    RemoteEndpointExecutor,
    _client_status_state,
)


class _Transport:
    host = "cluster"



def _run(native_id, *, backend="psij"):
    endpoint = endpoint_for_executor(RemoteEndpointExecutor(_Transport(), backend=backend))
    return RunRecord(
        id=native_id or "unconfirmed", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        task="DOS",
        execution=Execution(
            endpoint, WorkingDirectory("/local/case", "/work/case"), native_id,
        ),
        files=FileManifest(),
    )


def test_remote_endpoint_status_queries_opaque_handles_individually():
    executor = RemoteEndpointExecutor(_Transport(), backend="psij")
    calls = []
    replies = {
        "job-1": {"state": "running", "return_code": None},
        "job-2": {
            "state": "failed",
            "return_code": None,
            "failure_phase": "launch",
            "error": "FileNotFoundError('missing executable')",
        },
    }

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return replies[args[0]]

    executor._json_request = request
    runs = (_run("job-1"), _run("job-2"))

    statuses = executor.status(*runs)

    assert calls == [
        ("status", ("job-1",), None),
        ("status", ("job-2",), None),
    ]
    assert tuple(status.backend_id for status in statuses) == ("job-1", "job-2")
    assert tuple(status.state for status in statuses) == (
        RunState.RUNNING,
        RunState.FAILED,
    )
    assert statuses[1].return_code is None
    assert statuses[1].failure_phase == FailurePhase.LAUNCH
    assert "missing executable" in statuses[1].message


def test_direct_remote_status_uses_opaque_handle_without_client_receipt_path():
    executor = RemoteEndpointExecutor(_Transport(), backend="local")
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {"state": "finished", "return_code": 0}

    executor._json_request = request
    run = _run("opaque-process", backend="local")
    (status,) = executor.status(run)

    assert status.state == RunState.DOWNLOADING
    assert calls == [
        (
            "status",
            ("opaque-process",),
            None,
        )
    ]


def test_remote_endpoint_status_rejects_reference_for_other_backend():
    executor = RemoteEndpointExecutor(_Transport(), backend="local")

    with pytest.raises(ValueError, match="different remote endpoint"):
        executor.status(_run("job", backend="psij"))


def test_remote_recovery_passes_submission_timestamp_to_stateless_service():
    executor = RemoteEndpointExecutor(_Transport())
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {
            "state": "running",
            "remote_id": "opaque-job",
            "directory": "/work/case",
        }

    executor._json_request = request
    submission = _run(None, backend="local")

    status = executor.recover_submission(
        "run-1",
        submission,
        submitted_at="2026-10-13T08:00:00+02:00",
    )

    assert calls == [
        (
            "recover",
            ("/work/case", "run-1"),
            {
                "executor": "local",
                "executor_options": {},
                "submitted_at": "2026-10-13T08:00:00+02:00",
            },
        )
    ]
    assert status.backend_id == "opaque-job"
    assert status.state == RunState.RUNNING


def test_remote_recovery_can_report_unsupported_backend_recovery():
    executor = RemoteEndpointExecutor(_Transport())

    async def request(action, *args, payload=None):
        return {"unsupported": True, "error": "identity cannot be reconstructed"}

    executor._json_request = request

    with pytest.raises(NotImplementedError, match="cannot be reconstructed"):
        executor.recover_submission(
            "run-1", _run(None, backend="local")
        )


@pytest.mark.parametrize(
    "backend_state, phase, expected",
    [
        ("running", None, RunState.RUNNING),
        ("finished", None, RunState.DOWNLOADING),
        ("failed", "execution", RunState.DOWNLOADING),
        ("failed", "download", RunState.DOWNLOADING),
        ("failed", "launch", RunState.FAILED),
        ("failed", "submission", RunState.FAILED),
        ("cancelled", None, RunState.CANCELLED),
        ("unknown", None, RunState.UNKNOWN),
    ],
)
def test_backend_status_is_mapped_to_client_owned_lifecycle(
    backend_state, phase, expected
):
    state = {"state": backend_state, "failure_phase": phase}

    assert _client_status_state(state) == expected
