import pytest

import ase2sprkkr.execution.manager as manager_module
from ase2sprkkr.execution import (
    ArtifactAvailabilityState,
    ExecutionManager,
    FailurePhase,
    FileManifest,
    FileEntry,
    RunState,
    RunStore,
)


def test_result_snapshot_sync_failure_is_artifact_failure(tmp_path, monkeypatch):
    output = tmp_path / "result.out"
    output.write_text("result")
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        result_files=ResultFilesRecord(
            output_file=str(output),
            output_files=(str(output),),
            availability_state=ArtifactAvailabilityState.PENDING,
        ),
    )
    record = manager.transition(record, state=RunState.POSTPROCESSING)

    def fail_sync(_path):
        raise OSError("cannot synchronize result snapshot")

    monkeypatch.setattr(manager_module, "sync_file", fail_sync)

    failed = manager.complete(
        record,
        state=RunState.FINISHED,
        return_code=0,
        artifacts_available=True,
    )

    assert failed.state == RunState.FAILED
    assert failed.failure_phase == FailurePhase.ARTIFACTS
    assert failed.return_code == 0
    assert failed.finished_at is not None
    assert failed.files.captured_at is None
    assert failed.files.availability_state == ArtifactAvailabilityState.FAILED
    assert "cannot synchronize result snapshot" in failed.message
