from ase2sprkkr.execution import RunState, SubmissionKey
from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"

    def __init__(self, signature):
        self.signature = signature
        self.paths = []

    async def file_signature(self, path):
        self.paths.append(path)
        return self.signature


def _submission(directory="/scratch/case"):
    return SubmissionKey(
        endpoint=ExecutorEndpoint(executor="ssh", host="cluster"),
        directory=directory,
        task="DOS",
        dataset="case",
    )


def test_remote_recovery_returns_none_only_for_absent_helper_state():
    transport = _Transport(None)
    executor = RemoteEndpointExecutor(transport)

    status = executor.recover_submission("run-123", _submission())

    assert status is None
    assert transport.paths == [
        "/scratch/case/.ase2sprkkr/runs/run-123/state.json"
    ]


def test_remote_recovery_reconstructs_reference_and_status():
    transport = _Transport((10, 20))
    executor = RemoteEndpointExecutor(transport)

    async def status(directory, run_id):
        assert directory == "/scratch/case"
        assert run_id == "run-123"
        return {
            "state": "running",
            "return_code": None,
            "error": None,
            "input_file": "/scratch/case/run.inp",
            "output_file": "/scratch/case/run.out",
            "error_file": "/scratch/case/run.err",
        }

    executor._status = status

    recovered = executor.recover_submission("run-123", _submission())

    assert recovered.state == RunState.RUNNING
    assert recovered.reference.native_id == "run-123"
    assert recovered.reference.executor == "remote:local"
    assert recovered.reference.computer == "cluster"
    assert recovered.reference.directory == "/scratch/case"
    assert recovered.reference.output_file == "/scratch/case/run.out"
