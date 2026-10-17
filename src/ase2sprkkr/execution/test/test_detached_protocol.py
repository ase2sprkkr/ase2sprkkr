from io import BytesIO

from ase2sprkkr.execution import RunState
from ase2sprkkr.execution.detached import (
    _run_directory,
    _write_state,
    clean_detached_runs,
    follow_detached,
)


def test_public_execution_failed_state_is_terminal_for_follow_and_cleanup(tmp_path):
    run_id = "psij-failed"
    _run_directory(tmp_path, run_id).mkdir(parents=True)
    _write_state(
        tmp_path,
        run_id,
        id=run_id,
        state=RunState.EXECUTION_FAILED.value,
        output_file=None,
        error_file=None,
        finished_at="2026-10-16T12:00:00+00:00",
    )

    state = follow_detached(
        tmp_path,
        run_id,
        BytesIO(),
        BytesIO(),
        poll_interval=0,
    )

    assert state["state"] == RunState.EXECUTION_FAILED.value
    assert clean_detached_runs(tmp_path)["removed"] == [run_id]
