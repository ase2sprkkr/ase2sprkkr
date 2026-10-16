from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"


def _reference(native_id, *, backend="psij", context=None):
    return RunReference(
        executor=f"remote:{backend}",
        native_id=native_id,
        computer="cluster",
        directory="/work/case",
        context={} if context is None else context,
    )


def test_remote_endpoint_status_queries_opaque_ids():
    executor = RemoteEndpointExecutor(_Transport(), backend="psij")
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        remote_id = args[0]
        if remote_id == "job-1":
            return {"state": "running", "return_code": None}
        return {"state": "finished", "return_code": 0}

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
        RunState.FINISHED,
    )
    assert statuses[1].return_code == 0


def test_remote_endpoint_status_propagates_backend_error():
    executor = RemoteEndpointExecutor(_Transport(), backend="psij")

    async def request(action, *args, payload=None):
        raise OSError("scheduler unavailable")

    executor._json_request = request

    with pytest.raises(OSError, match="scheduler unavailable"):
        executor.status(_reference("job-1"))


def test_remote_recovery_accepts_submission_timestamp_and_context(monkeypatch):
    executor = RemoteEndpointExecutor(_Transport())
    calls = []

    async def recover(run_id, submission, *, submitted_at=None):
        calls.append((run_id, submission.directory, submitted_at, submission.context))
        return None

    monkeypatch.setattr(executor, "_recover_submission", recover)
    submission = SimpleNamespace(
        directory="/work/case",
        context={"launch": {"remote_root": "/work"}},
    )

    assert executor.recover_submission(
        "run-1",
        submission,
        submitted_at="2026-10-13T08:00:00+02:00",
    ) is None
    assert calls == [
        (
            "run-1",
            "/work/case",
            "2026-10-13T08:00:00+02:00",
            {"launch": {"remote_root": "/work"}},
        )
    ]
