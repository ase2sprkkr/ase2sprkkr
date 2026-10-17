from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import RunReference, RunState, RunStatus
from ase2sprkkr.execution import remote_service


def test_direct_remote_process_cannot_recover_lost_submit_identity():
    with pytest.raises(NotImplementedError, match="cannot recover"):
        remote_service.recover_remote(
            "client-run",
            {"executor": "local", "directory": "/scratch/case"},
        )


def test_psij_recovery_uses_backend_native_identity(monkeypatch):
    calls = []

    class Backend:
        backend = "slurm"
        url = None

        def recover_submission(self, run_id, submission, *, submitted_at=None):
            calls.append((run_id, submission.directory, submitted_at))
            reference = RunReference(
                executor="psij:slurm",
                native_id="scheduler-42",
                directory="/scratch/case",
                input_file="/scratch/case/run.inp",
                output_file="/scratch/case/run.out",
                error_file="/scratch/case/run.err",
            )
            return RunStatus(reference, RunState.RUNNING)

    monkeypatch.setattr(remote_service, "_psij_executor", lambda options: Backend())

    recovered = remote_service.recover_remote(
        "client-run",
        {
            "executor": "psij",
            "executor_options": {"backend": "slurm", "url": None},
            "directory": "/scratch/case",
        },
        submitted_at="2026-10-13T08:00:00+02:00",
    )

    assert calls == [
        ("client-run", "/scratch/case", "2026-10-13T08:00:00+02:00")
    ]
    assert recovered["state"] == RunState.RUNNING.value
    assert recovered["directory"] == "/scratch/case"
    assert recovered["output_file"] == "/scratch/case/run.out"

    handle = remote_service._decode_handle(recovered["remote_id"])
    assert handle == {
        "kind": "psij",
        "backend": "slurm",
        "url": None,
        "native_id": "scheduler-42",
    }


def test_psij_recovery_returns_none_when_backend_has_no_match(monkeypatch):
    class Backend:
        backend = "slurm"
        url = None

        def recover_submission(self, run_id, submission, *, submitted_at=None):
            return None

    monkeypatch.setattr(remote_service, "_psij_executor", lambda options: Backend())

    assert remote_service.recover_remote(
        "client-run",
        {
            "executor": "psij",
            "executor_options": {"backend": "slurm"},
            "directory": "/scratch/case",
        },
    ) is None
