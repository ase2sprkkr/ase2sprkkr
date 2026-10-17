import asyncio
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import ExecutionSpec, FailurePhase, RunReference, RunState
from ase2sprkkr.execution.base import SubmissionNotStarted
from ase2sprkkr.execution.remote import RemoteRunner
from ase2sprkkr.execution.remote_endpoint import _RemoteRun


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None)


class _Transport:
    def __init__(self):
        self.signature_calls = []

    async def run(self, argv):
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    async def file_signature(self, path):
        self.signature_calls.append(path)
        return None


class _Executor:
    backend = "local"
    computer = "cluster"

    def __init__(self):
        self.transport = _Transport()

    async def _status(self, directory, run_id):
        return {
            "state": "cancelled",
            "started_at": "execution-time",
        }


def _run(tmp_path, executor=None):
    executor = executor or _Executor()
    return _RemoteRun(
        executor,
        _Reader(),
        ExecutionSpec(
            argv=(),
            input_file=None,
            output_file=str(tmp_path / "calculation.out"),
            directory=str(tmp_path),
        ),
        reference=RunReference(
            executor="remote:local",
            native_id="job-cancelled",
            computer="cluster",
            directory="/work/case",
            output_file="/work/case/calculation.out",
        ),
        state=RunState.CANCELLED,
    )


def test_remote_post_output_handling_is_allowed_after_running(tmp_path, monkeypatch):
    executor = _Executor()
    run = _run(tmp_path, executor)

    async def base_monitor(self):
        return self.result

    monkeypatch.setattr(RemoteRunner, "_monitor", base_monitor)
    asyncio.run(run._monitor())

    assert executor.transport.signature_calls == ["/work/case/calculation.out"]


def test_remote_helper_failure_preserves_launch_phase(tmp_path, monkeypatch):
    executor = _Executor()
    run = _run(tmp_path, executor)
    run._reference = None
    run._run_id = "remote-launch-failure"

    async def map_directory(directory):
        return "/work/case"

    async def prepare_request(run_id, remote_directory):
        return {"id": run_id, "directory": remote_directory}

    async def fail_submit(request):
        raise SubmissionNotStarted(
            "Remote detached worker launch failed: OSError('missing helper')"
        )

    executor.map_directory = map_directory
    executor._submit = fail_submit
    monkeypatch.setattr(run, "_prepare_request", prepare_request)

    with pytest.raises(SubmissionNotStarted, match="missing helper"):
        asyncio.run(run._submit_async())

    assert run.failure_phase == FailurePhase.LAUNCH
    assert "missing helper" in run.failure_message
