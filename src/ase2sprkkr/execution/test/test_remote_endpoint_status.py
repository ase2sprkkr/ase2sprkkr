import asyncio
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import ExecutionSpec, FailurePhase, RunReference, RunState
from ase2sprkkr.execution.base import SubmissionNotStarted
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor, _RemoteRun


class _Transport:
    host = "cluster"


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None)


def test_remote_endpoint_status_uses_pooled_endpoint_backend():
    executor = RemoteEndpointExecutor(_Transport(), backend="psij")
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return [
            {"state": "running", "return_code": None},
            {
                "state": "failed",
                "return_code": None,
                "failure_phase": "launch",
                "error": "FileNotFoundError('missing executable')",
            },
        ]

    executor._json_request = request
    references = (
        RunReference(
            executor="remote:psij",
            native_id="job-1",
            computer="cluster",
            directory="/work/case",
        ),
        RunReference(
            executor="remote:psij",
            native_id="job-2",
            computer="cluster",
            directory="/work/case",
        ),
    )

    statuses = executor.status(*references)

    assert calls == [("status", ("/work/case", "job-1", "job-2"), None)]
    assert tuple(status.reference for status in statuses) == references
    assert tuple(status.state for status in statuses) == (
        RunState.RUNNING,
        RunState.EXECUTION_FAILED,
    )
    assert statuses[1].return_code is None
    assert statuses[1].failure_phase == FailurePhase.LAUNCH
    assert "missing executable" in statuses[1].message


def test_remote_endpoint_rejects_missing_status_rows():
    executor = RemoteEndpointExecutor(_Transport())

    async def request(action, *args, payload=None):
        return [{"state": "running"}]

    executor._json_request = request
    references = tuple(
        RunReference(
            executor="remote:local",
            native_id=f"job-{index}",
            computer="cluster",
            directory="/work/case",
        )
        for index in range(2)
    )

    with pytest.raises(RuntimeError, match="wrong number"):
        executor.status(*references)


def test_remote_recovery_accepts_submission_timestamp(monkeypatch):
    executor = RemoteEndpointExecutor(_Transport())
    calls = []

    async def recover(run_id, submission):
        calls.append((run_id, submission.directory))
        return None

    monkeypatch.setattr(executor, "_recover_submission", recover)

    assert executor.recover_submission(
        "run-1",
        SimpleNamespace(directory="/work/case"),
        submitted_at="2026-10-13T08:00:00+02:00",
    ) is None
    assert calls == [("run-1", "/work/case")]


def test_remote_detached_helper_failure_includes_remote_reason():
    executor = RemoteEndpointExecutor(_Transport())

    async def request(action, *args, payload=None):
        assert action == "submit"
        return {
            "state": "failed",
            "failure_phase": "launch",
            "error": "OSError('remote helper missing')",
        }

    executor._json_request = request

    with pytest.raises(
        SubmissionNotStarted,
        match="Remote detached worker launch failed.*remote helper missing",
    ):
        asyncio.run(executor._submit({"directory": "/work/case"}))


def _remote_run_with_statuses(tmp_path, statuses):
    class Executor:
        async def _status(self, directory, run_id):
            return statuses.pop(0)

    reference = RunReference(
        executor="remote:psij",
        native_id="job-1",
        computer="cluster",
        directory="/work/case",
    )
    run = _RemoteRun(
        Executor(),
        _Reader(),
        ExecutionSpec(
            argv=(),
            input_file=None,
            output_file=None,
            directory=str(tmp_path),
        ),
        reference=reference,
        state=RunState.QUEUED,
    )
    return run


def test_remote_waits_for_running_evidence(tmp_path, monkeypatch):
    statuses = [
        {"state": "queued", "submitted_at": "submission-time", "started_at": None},
        {"state": "running", "started_at": "execution-time"},
    ]
    run = _remote_run_with_statuses(tmp_path, statuses)

    async def no_sleep(_delay):
        return None

    monkeypatch.setattr("ase2sprkkr.execution.remote.asyncio.sleep", no_sleep)
    status = asyncio.run(run._wait_until_started())

    assert status["state"] == "running"
    assert status["started_at"] == "execution-time"
    assert statuses == []


def test_remote_cancel_before_running_is_not_treated_as_started(tmp_path):
    run = _remote_run_with_statuses(
        tmp_path,
        [{"state": "cancelled", "submitted_at": "submission-time", "started_at": None}],
    )

    status = asyncio.run(run._wait_until_started())

    assert status["state"] == "cancelled"
    assert status["started_at"] is None


def test_remote_cancel_after_running_keeps_execution_evidence(tmp_path):
    run = _remote_run_with_statuses(
        tmp_path,
        [{"state": "cancelled", "started_at": "execution-time"}],
    )

    status = asyncio.run(run._wait_until_started())

    assert status["state"] == "cancelled"
    assert status["started_at"] == "execution-time"


def test_remote_cancel_before_running_skips_post_output_handling(tmp_path):
    class Transport:
        async def file_signature(self, path):
            raise AssertionError("output must not be inspected before RUNNING")

    class Executor:
        transport = Transport()

        async def _status(self, directory, run_id):
            return {
                "state": "cancelled",
                "submitted_at": "submission-time",
                "started_at": None,
            }

    run = _RemoteRun(
        Executor(),
        _Reader(),
        ExecutionSpec(
            argv=(),
            input_file=None,
            output_file=str(tmp_path / "calculation.out"),
            directory=str(tmp_path),
        ),
        reference=RunReference(
            executor="remote:psij",
            native_id="job-cancelled",
            computer="cluster",
            directory="/work/case",
            output_file="/work/case/calculation.out",
        ),
        state=RunState.QUEUED,
    )

    result = asyncio.run(run._monitor())

    assert result.return_code is None
    assert run.state == RunState.CANCELLED
    assert not (tmp_path / "calculation.out").exists()
