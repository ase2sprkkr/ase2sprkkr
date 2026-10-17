from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import FailurePhase, RunReference, RunState
from ase2sprkkr.execution.remote_endpoint import (
    RemoteEndpointExecutor,
    _client_status_state,
)


class _Transport:
    host = "cluster"



def _reference(native_id, *, backend="psij"):
    return RunReference(
        executor=f"remote:{backend}",
        native_id=native_id,
        computer="cluster",
        directory="/work/case",
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
    references = (_reference("job-1"), _reference("job-2"))

    statuses = executor.status(*references)

    assert calls == [
        ("status", ("job-1",), None),
        ("status", ("job-2",), None),
    ]
    assert tuple(status.reference for status in statuses) == references
    assert tuple(status.state for status in statuses) == (
        RunState.RUNNING,
        RunState.FAILED,
    )
    assert statuses[1].return_code is None
    assert statuses[1].failure_phase == FailurePhase.LAUNCH
    assert "missing executable" in statuses[1].message


def test_direct_remote_status_passes_client_owned_receipt_location():
    executor = RemoteEndpointExecutor(_Transport(), backend="local")
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {"state": "finished", "return_code": 0}

    executor._json_request = request
    reference = RunReference(
        executor="remote:local",
        native_id="opaque-process",
        directory="/work/case",
        context={},
    )
    (status,) = executor.status(reference)

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

    with pytest.raises(ValueError, match="belongs to"):
        executor.status(_reference("job", backend="psij"))


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
    submission = SimpleNamespace(directory="/work/case", context=None)

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
    assert status.reference.native_id == "opaque-job"
    assert status.state == RunState.RUNNING


def test_remote_recovery_can_report_unsupported_backend_recovery():
    executor = RemoteEndpointExecutor(_Transport())

    async def request(action, *args, payload=None):
        return {"unsupported": True, "error": "identity cannot be reconstructed"}

    executor._json_request = request

    with pytest.raises(NotImplementedError, match="cannot be reconstructed"):
        executor.recover_submission(
            "run-1", SimpleNamespace(directory="/work/case", context=None)
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
