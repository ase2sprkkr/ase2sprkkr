from types import SimpleNamespace

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"



def test_remote_endpoint_status_uses_pooled_endpoint_backend():
    executor = RemoteEndpointExecutor(_Transport(), backend="psij")
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return [
            {"state": "running", "return_code": None},
            {"state": "finished", "return_code": 0},
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
        RunState.FINISHED,
    )
    assert statuses[1].return_code == 0


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
