import sys
import time

from ase2sprkkr.execution import RunState
from ase2sprkkr.execution.remote_service import cancel_remote, status_remote, submit_remote


def _submitted(tmp_path, code):
    output = tmp_path / "out"
    return submit_remote(
        {
            "id": "cancel-run",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [sys.executable, "-c", code],
            "input_file": None,
            "stdin": None,
            "output_file": str(output),
            "error_file": str(output) + ".err",
            "launch": {},
            "mpi_argv": None,
        }
    )


def _finished(remote_id):
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        state = status_remote(remote_id)
        if state["state"] != RunState.RUNNING.value:
            return state
        time.sleep(0.02)
    raise AssertionError("remote process did not finish")


def test_late_cancel_preserves_completed_process_state(tmp_path):
    submitted = _submitted(tmp_path, "print('done')")
    before = _finished(submitted["remote_id"])

    after = cancel_remote(submitted["remote_id"])

    assert before["state"] == RunState.FINISHED.value
    assert before["return_code"] == 0
    assert after["state"] == RunState.FINISHED.value
    assert after["return_code"] == 0


def test_cancel_live_process_reports_cancelled(tmp_path):
    submitted = _submitted(tmp_path, "import time; time.sleep(30)")

    state = cancel_remote(submitted["remote_id"], timeout=0.2)

    assert state["state"] == RunState.CANCELLED.value
