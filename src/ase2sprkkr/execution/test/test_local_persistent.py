"""Local persistent execution has no independent lifecycle database."""

import os
import sys
import time
from dataclasses import replace

import pytest

from ase2sprkkr.execution import (
    Execution, ExecutionSpec, ExecutorEndpoint, FileEntry, FileManifest,
    LocalExecutor, RunRecord, RunState, WorkingDirectory,
)
from ase2sprkkr.execution import local_backend


pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="reconnectable local identity uses Linux /proc start ticks",
)


def _spec(tmp_path, command=None):
    return ExecutionSpec(
        argv=command or (sys.executable, "-c", "print('output', flush=True)"),
        input_file=None, output_file=str(tmp_path / "out.txt"),
        error_file=str(tmp_path / "err.txt"), directory=str(tmp_path),
    )


def _record(tmp_path, backend_id):
    return RunRecord(
        id="run", state=RunState.RUNNING, task="DOS",
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(str(tmp_path)), backend_id,
        ),
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "out.txt")),
            stderr=FileEntry(local_path=str(tmp_path / "err.txt")),
        ),
    )


def test_run_reconnect_and_immutable_backend_facts(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_persistent_run(_spec(tmp_path), run_id="run")
    runner.submit()
    identity = runner.backend_id
    assert identity.startswith("local:run:")
    assert runner.wait().return_code == 0
    assert (tmp_path / "out.txt").read_text() == "output\n"
    assert executor.status(_record(tmp_path, identity))[0].state == RunState.FINISHED
    assert executor.recover_submission("run", _record(tmp_path, None)).backend_id == identity
    metadata = tmp_path / ".ase2sprkkr" / "runs" / "run"
    assert sorted(p.name for p in metadata.iterdir()) == ["identity.json", "outcome.json"]

    seen = []
    attached = executor.attach(
        _record(tmp_path, identity),
        output_listener=lambda text, kind: seen.append((text, kind)),
    )
    assert attached.wait().return_code == 0
    assert seen == [("output\n", "out")]


def test_failed_process_preserves_execution_failure(tmp_path):
    runner = LocalExecutor().create_persistent_run(
        _spec(tmp_path, (sys.executable, "-c", "import sys; sys.exit(7)")),
        run_id="fail",
    )
    runner.submit()
    assert runner.wait().return_code == 7
    assert runner.state == RunState.FAILED
    assert runner.failure_phase.value == "execution"


def test_reused_pid_is_not_a_running_process(tmp_path):
    runner = LocalExecutor().create_persistent_run(_spec(tmp_path), run_id="pid")
    runner.submit()
    runner.wait()
    identity = local_backend._parse_id(runner.backend_id)
    assert local_backend._alive({**identity, "start": "0"}) is False


def test_cancel_does_not_block_and_publishes_receipt(tmp_path):
    runner = LocalExecutor().create_persistent_run(
        _spec(tmp_path, (
            sys.executable, "-c",
            "import time; print('started', flush=True); time.sleep(30)",
        )),
        run_id="cancel",
    )
    runner.submit()
    started = time.monotonic()
    runner.cancel()
    assert time.monotonic() - started < 1
    runner.wait()
    assert runner.state == RunState.CANCELLED
    assert local_backend.query(tmp_path, runner.backend_id).state == RunState.CANCELLED


def test_missing_worker_outcome_is_unknown_not_success(tmp_path):
    record = _record(tmp_path, "local:missing:9999999:1:boot")
    status = LocalExecutor().status(record)[0]
    assert status.state == RunState.UNKNOWN


def test_wildcard_files_have_individual_status(tmp_path):
    from ase2sprkkr.execution.run_model import FileManifest
    (tmp_path / "one.dat").write_text("x")
    manifest = FileManifest(outputs=(
        FileEntry(local_path=str(tmp_path / "*.dat"), availability_state="pending"),
    ))
    identity = local_backend._identity("wildcard")
    from ase2sprkkr.execution._durability import durable_mkdir, durable_write_json
    metadata = local_backend._run_dir(tmp_path, "wildcard")
    durable_mkdir(metadata, parents=True)
    durable_write_json(metadata / "identity.json", identity)
    durable_write_json(metadata / "outcome.json", {
        "return_code": 0, "phase": None, "message": None,
        "cancelled": False, "published": [True, True],
    })
    updated = local_backend.local_manifest(manifest, tmp_path, local_backend._backend_id(identity))
    assert len(updated.outputs) == 2
    assert updated.outputs[0].availability_state == "available"
    assert updated.outputs[1].source_pattern == str(tmp_path / "*.dat")
    assert updated.outputs[1].availability_state == "available"
