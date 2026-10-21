"""Execution outcome and file availability are independent persisted facts."""

import pytest

from ase2sprkkr.execution import (
    ArtifactAvailabilityState, FailurePhase, FileEntry,
    FileManifest, RunState, RunStore,
)

from ase2sprkkr.execution.manager import ExecutionManager


def _record(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "run.out"), availability_state="pending"),
        ),
    )
    record = manager.transition(record, state=RunState.POSTPROCESSING)
    return manager, record


@pytest.mark.parametrize(
    "requested_state, return_code, phase, backend_error, transfer_error, "
    "artifacts_available, expected_state, expected_phase",
    [
        pytest.param(
            RunState.FINISHED, 0, None, None, "SFTP connection lost",
            True, RunState.FAILED, FailurePhase.DOWNLOAD,
            id="successful-backend-failed-transfer",
        ),
        pytest.param(
            RunState.FAILED, 17, FailurePhase.EXECUTION,
            "SPR-KKR returned 17", "diagnostic files unavailable",
            False, RunState.FAILED, FailurePhase.EXECUTION,
            id="execution-failure-overrides-transfer-error",
        ),
        pytest.param(
            RunState.CANCELLED, None, None, None, "download interrupted",
            False, RunState.CANCELLED, None,
            id="cancellation-preserves-outcome",
        ),
    ],
)
def test_backend_outcome_and_transfer_availability_are_independent(
    tmp_path, requested_state, return_code, phase, backend_error,
    transfer_error, artifacts_available, expected_state, expected_phase,
):
    manager, record = _record(tmp_path)
    failed_files = record.files.with_file_availability(
        str(tmp_path / "run.out"), "failed", transfer_error,
    )
    saved = manager.complete(
        record, state=requested_state, return_code=return_code,
        failure_phase=phase, message=backend_error,
        artifacts_available=artifacts_available, files=failed_files,
    )

    assert manager.store.get(record.id) == saved
    assert saved.state is expected_state
    assert saved.failure_phase is expected_phase
    assert saved.return_code == return_code
    assert saved.message == (
        transfer_error if expected_phase is FailurePhase.DOWNLOAD
        else backend_error
    )
    assert saved.files.availability_state == ArtifactAvailabilityState.FAILED
    assert saved.files.availability_message == transfer_error
    assert saved.files.captured_at is None


def test_successful_run_captures_available_files(tmp_path):
    manager, record = _record(tmp_path)
    (tmp_path / "run.out").write_text("complete")
    available = record.files.with_file_availability(str(tmp_path / "run.out"), "available")
    saved = manager.complete(
        record, state=RunState.FINISHED, return_code=0,
        artifacts_available=True, files=available,
    )
    assert saved.state is RunState.FINISHED
    assert saved.return_code == 0
    assert saved.files.availability_state == ArtifactAvailabilityState.AVAILABLE.value
    assert saved.files.captured_at is not None
    assert len(saved.files.output_versions) == 1


def test_resume_download_rejects_execution_failure_and_nonzero_code(tmp_path):
    manager, record = _record(tmp_path)
    bad = manager.complete(
        record, state=RunState.FAILED, return_code=17,
        failure_phase=FailurePhase.EXECUTION, artifacts_available=False,
    )
    with pytest.raises(ValueError, match="Only failed downloads"):
        manager.resume_download(bad.id)


def test_retry_failed_download_preserves_backend_then_captures_clean_snapshot(
    tmp_path,
):
    manager, initial = _record(tmp_path)
    output = tmp_path / "run.out"
    output.write_bytes(b"partial")
    initial = manager.transition(
        initial, state=RunState.DOWNLOADING,
        execution=initial.execution.submitted("existing-backend"),
    )
    failed_files = initial.files.with_file_availability(
        str(output), "failed", "interrupted transfer",
    )
    failed = manager.complete(
        initial, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.DOWNLOAD,
        message="interrupted transfer",
        files=failed_files, artifacts_available=False,
    )

    resumed = manager.resume_download(initial.id)

    assert resumed.state is RunState.DOWNLOADING
    assert resumed.execution == failed.execution
    assert resumed.execution.backend_id == "existing-backend"
    assert resumed.files.stdout.availability_state == "pending"
    assert resumed.files.stdout.availability_message is None
    assert resumed.finished_at is None
    assert resumed.failure_phase is None
    assert output.read_bytes() == b"partial"
    assert manager.store.get(initial.id) == resumed

    output.write_bytes(b"complete")
    available = resumed.files.with_file_availability(str(output), "available")
    complete = manager.complete(
        resumed, state=RunState.FINISHED, return_code=0,
        files=available, artifacts_available=True,
    )
    assert complete.state is RunState.FINISHED
    assert complete.files.failed_artifacts == ()
    assert complete.files.captured_at is not None
    assert len(complete.files.output_versions) == 1
    assert manager.output_file_is_current(complete)


def test_resume_download_retains_previously_available_artifacts(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "stdout"
    other = tmp_path / "other"
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(
            stdout=FileEntry(local_path=str(output), availability_state="available"),
            outputs=(FileEntry(local_path=str(other), availability_state="failed",
                               availability_message="temporary transport error"),),
        ),
    )
    record = manager.transition(record, state=RunState.DOWNLOADING,
                                execution=record.execution.submitted("existing-backend"))
    manager.complete(record, state=RunState.FAILED, return_code=0,
                     failure_phase=FailurePhase.DOWNLOAD, artifacts_available=False)
    retried = manager.resume_download(record.id)
    assert retried.files.stdout.availability_state == "available"
    assert retried.files.outputs[0].availability_state == "pending"
    assert retried.files.outputs[0].availability_message is None
