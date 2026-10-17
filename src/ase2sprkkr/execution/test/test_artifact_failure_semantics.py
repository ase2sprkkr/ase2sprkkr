"""Execution outcome and file availability are independent persisted facts."""

from dataclasses import replace

from ase2sprkkr.execution import (
    ArtifactAvailabilityState, ExecutionManager, FailurePhase, FileEntry,
    FileManifest, RunState, RunStore,
)


def _record(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "run.out")),
            availability_state="pending",
        ),
    )
    record = manager.transition(record, state=RunState.POSTPROCESSING)
    return manager, record


def test_successful_backend_with_failed_transfer_is_failed_run(tmp_path):
    manager, record = _record(tmp_path)
    failed = replace(
        record.files,
        availability_state="failed",
        availability_message="SFTP connection lost",
    )
    result = manager.complete(
        record,
        state=RunState.FINISHED,
        return_code=0,
        artifacts_available=True,
        files=failed,
    )
    saved = manager.store.get(record.id)
    assert result == saved
    assert saved.state is RunState.FAILED
    assert saved.failure_phase is FailurePhase.DOWNLOAD
    assert saved.return_code == 0
    assert saved.message == "SFTP connection lost"
    assert saved.files.availability_state == ArtifactAvailabilityState.FAILED
    assert saved.files.captured_at is None


def test_execution_failure_keeps_phase_and_message_when_transfer_fails(tmp_path):
    manager, record = _record(tmp_path)
    files = replace(
        record.files,
        availability_state="failed",
        availability_message="diagnostic files unavailable",
    )
    saved = manager.complete(
        record,
        state=RunState.FAILED,
        failure_phase=FailurePhase.EXECUTION,
        message="SPR-KKR returned 17",
        return_code=17,
        artifacts_available=False,
        files=files,
    )
    assert saved.state is RunState.FAILED
    assert saved.failure_phase is FailurePhase.EXECUTION
    assert saved.message == "SPR-KKR returned 17"
    assert saved.return_code == 17
    assert saved.files.availability_state == ArtifactAvailabilityState.FAILED
    assert saved.files.availability_message == "diagnostic files unavailable"


def test_successful_run_captures_available_files(tmp_path):
    manager, record = _record(tmp_path)
    (tmp_path / "run.out").write_text("complete")
    available = replace(record.files, availability_state="available")
    saved = manager.complete(
        record, state=RunState.FINISHED, return_code=0,
        artifacts_available=True, files=available,
    )
    assert saved.state is RunState.FINISHED
    assert saved.return_code == 0
    assert saved.files.availability_state == ArtifactAvailabilityState.AVAILABLE
    assert saved.files.captured_at is not None
    assert len(saved.files.output_versions) == 1


def test_cancellation_preserves_outcome_when_diagnostic_files_fail(tmp_path):
    manager, record = _record(tmp_path)
    files = replace(
        record.files, availability_state="failed",
        availability_message="download interrupted",
    )
    saved = manager.complete(
        record, state=RunState.CANCELLED, files=files,
        artifacts_available=False,
    )
    assert saved.state is RunState.CANCELLED
    assert saved.failure_phase is None
    assert saved.files.availability_state == ArtifactAvailabilityState.FAILED
