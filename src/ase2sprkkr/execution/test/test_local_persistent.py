import os
import sys

import pytest

from ase2sprkkr.execution import (
    ExecutionSpec, LocalExecutor, RunState, RunRecord,
    Execution, ExecutorEndpoint, WorkingDirectory, FileManifest, FileEntry,
)


pytestmark = pytest.mark.skipif(
    os.name == "nt",
    reason="the detached helper currently uses POSIX file locking",
)


def _spec(tmp_path):
    return ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import time; print('before', flush=True); "
            "time.sleep(0.1); print('after', flush=True)",
        ),
        input_file=None,
        output_file=str(tmp_path / "calculation.out"),
        directory=str(tmp_path),
    )


def _record(tmp_path, backend_id):
    return RunRecord(
        id=backend_id,
        state=RunState.RUNNING,
        task="DOS",
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(str(tmp_path)), backend_id,
        ),
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "calculation.out")),
        ),
    )


def test_persistent_local_run_can_be_reattached_and_replayed(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_persistent_run(
        _spec(tmp_path), run_id="local-persistent",
    )
    runner.submit()
    assert runner.backend_id == "local-persistent"

    result = runner.wait()
    assert result.return_code == 0
    assert (tmp_path / "calculation.out").read_bytes() == b"before\nafter\n"
    record = _record(tmp_path, runner.backend_id)
    assert executor.status(record)[0].state == RunState.FINISHED

    # Reattach with only persisted execution and manifest information.
    (tmp_path / "calculation.out").write_text("stale\n")
    replay = []
    attached = LocalExecutor().attach(
        record,
        output_listener=lambda text, kind: replay.append((text, kind)),
        start=False,
    )
    assert replay == []
    attached.start()
    assert attached.wait().return_code == 0
    assert (tmp_path / "calculation.out").read_bytes() == b"before\nafter\n"
    assert "".join(text for text, kind in replay if kind == "out") == "before\nafter\n"


def test_persistent_local_submission_can_be_recovered_by_run_id(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_persistent_run(
        _spec(tmp_path), run_id="recover-local",
    )
    runner.submit()
    recovered = LocalExecutor().recover_submission(
        "recover-local", _record(tmp_path, None),
    )
    assert recovered is not None
    assert recovered.backend_id == "recover-local"
    assert recovered.state in {
        RunState.SUBMITTING, RunState.RUNNING, RunState.FINISHED,
    }
    runner.wait()
