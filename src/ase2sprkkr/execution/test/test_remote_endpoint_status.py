from dataclasses import replace

import pytest

from ase2sprkkr.execution import (
    Execution, ExecutorEndpoint, FileEntry, FileManifest, RunRecord,
    FailurePhase, RunState, WorkingDirectory,
)
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


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
    "facts, expected_state",
    [
        ({"state": "running"}, RunState.RUNNING),
        ({"state": "finished", "return_code": 0}, RunState.DOWNLOADING),
        ({"state": "failed", "return_code": 7, "failure_phase": "execution", "error": "job failed"}, RunState.DOWNLOADING),
        ({"state": "failed", "failure_phase": "download", "error": "transfer failed"}, RunState.DOWNLOADING),
        ({"state": "failed", "failure_phase": "launch", "error": "program missing"}, RunState.FAILED),
        ({"state": "failed", "failure_phase": "submission"}, RunState.FAILED),
        ({"state": "cancelled", "return_code": 130}, RunState.CANCELLED),
        ({"state": "unknown"}, RunState.UNKNOWN),
    ],
)
def test_remote_status_and_recovery_share_backend_facts_mapping(
    facts, expected_state,
):
    executor = RemoteEndpointExecutor(_Transport(), backend="local")
    calls = []

    async def request(action, *args, payload=None):
        calls.append(action)
        if action == "status":
            assert args == ("opaque-job",)
            return dict(facts)
        if action == "recover":
            assert args == ("/work/case", "run-id")
            return {
                **facts,
                "remote_id": "opaque-job",
                "directory": "/other/work",
                "output_file": "/other/work/run.out",
                "error_file": "/other/work/run.err",
            }
        pytest.fail(f"Unexpected operation: {action}")

    executor._json_request = request
    (queried,) = executor.status(_run("opaque-job", backend="local"))
    original = _run(None, backend="local")
    submission = replace(original, files=FileManifest(
        stdout=FileEntry(local_path="/local/case/run.out"),
        stderr=FileEntry(local_path="/local/case/run.err"),
    ))
    recovered = executor.recover_submission("run-id", submission)

    assert calls == ["status", "recover"]
    assert queried.state == recovered.state == expected_state
    assert queried.backend_id == recovered.backend_id == "opaque-job"
    assert queried.return_code == recovered.return_code == facts.get("return_code")
    assert queried.message == recovered.message == facts.get("error")
    expected_phase = facts.get("failure_phase")
    expected_phase = None if expected_phase is None else FailurePhase(expected_phase)
    assert queried.failure_phase == recovered.failure_phase == expected_phase
    assert queried.execution is None
    assert queried.files is None
    assert recovered.execution.directory.local == "/local/case"
    assert recovered.execution.directory.remote == "/other/work"
    assert recovered.files.stdout.local_path == "/local/case/run.out"
    assert recovered.files.stdout.remote_path == "/other/work/run.out"
    assert recovered.files.stderr.local_path == "/local/case/run.err"
    assert recovered.files.stderr.remote_path == "/other/work/run.err"


@pytest.mark.parametrize(
    "facts, expected_phase",
    [
        ({"return_code": 7, "failure_phase": "execution", "error": "failed"},
         FailurePhase.EXECUTION),
        ({"return_code": 0, "failure_phase": None, "error": None}, None),
        ({}, None),
    ],
)
def test_remote_runner_backend_facts_are_consistent(
    tmp_path, facts, expected_phase,
):
    """Backend metadata is applied to both the runner and its result."""
    from ase2sprkkr.execution.base import ExecutionSpec
    from ase2sprkkr.execution.remote import RemoteRunner

    runner = RemoteRunner(
        None, ExecutionSpec(
            argv=(), input_file=None, output_file=None,
            directory=str(tmp_path),
        ),
    )
    runner._apply_backend_facts(facts)

    assert runner.return_code == runner.result.return_code == facts.get("return_code")
    assert runner.failure_phase == expected_phase
    assert runner.failure_message == facts.get("error")
