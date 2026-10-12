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
