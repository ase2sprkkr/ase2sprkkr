import pytest

from ase2sprkkr.execution import RunState
from ase2sprkkr.execution import remote_service


def test_direct_remote_process_without_pending_identity_returns_none(tmp_path, monkeypatch):
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))
    assert remote_service.recover_remote(
        "a" * 32, {"executor": "local", "directory": "/scratch/case"},
    ) is None


def test_psij_recovery_uses_backend_native_identity(monkeypatch):
    calls = []

    class Backend:
        backend = "slurm"
        url = None

        def recover_native(self, run_id, directory, *, submitted_at=None):
            calls.append((run_id, directory, submitted_at))
            return {
                "native_id": "scheduler-42",
                "state": RunState.RUNNING,
                "directory": "/scratch/case",
                "stdout": "/scratch/case/run.out",
                "stderr": "/scratch/case/run.err",
                "return_code": None,
            }

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

        def recover_native(self, run_id, directory, *, submitted_at=None):
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
