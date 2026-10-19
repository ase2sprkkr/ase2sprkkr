from dataclasses import asdict

import pytest

from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.run_model import (
    CacheEntry, Execution, FileEntry, FileManifest, WorkingDirectory,
)


def test_execution_roundtrip_includes_backend_and_both_directories():
    original = Execution(
        ExecutorEndpoint(executor="ssh", host="cluster"),
        WorkingDirectory("/local/case", "/scratch/case"),
        backend_id="slurm-42",
    )
    assert Execution(**asdict(original)) == original
    assert original.directory.target == "/scratch/case"


def test_local_working_directory_does_not_need_remote_copy():
    run = Execution(ExecutorEndpoint(), WorkingDirectory("/tmp/work"))
    assert run.backend_id is None
    assert run.directory.target == "/tmp/work"


def test_manifest_roundtrip_includes_streams_and_cache_baseline():
    original = FileManifest(
        inputs=(FileEntry("/local/input", "/remote/input"),),
        outputs=(FileEntry("/local/result", "/remote/result"),),
        caches=(CacheEntry(FileEntry("/local/cache", "/remote/cache"), (100, 123)),),
        stdout=FileEntry("/local/out", "/remote/out"),
        stderr=FileEntry("/local/err", "/remote/err"),
    )
    restored = FileManifest(**asdict(original))
    assert restored == original
    assert restored.caches[0].before == (100, 123)


def test_empty_backend_identity_and_pathless_file_are_invalid():
    with pytest.raises(ValueError, match="backend_id"):
        Execution(ExecutorEndpoint(), WorkingDirectory("/tmp"), backend_id="")
    with pytest.raises(ValueError, match="location"):
        FileEntry()


def test_manifest_versions_and_availability_roundtrip():
    manifest = FileManifest(
        inputs=(FileEntry("/local/in", "/remote/in"),),
        outputs=(FileEntry("/local/data", "/remote/data"),),
        stdout=FileEntry("/local/stdout", "/remote/stdout"),
        stderr=FileEntry("/local/stderr", "/remote/stderr"),
    ).initialize_availability()
    restored = FileManifest(**asdict(manifest))
    assert restored == manifest
    assert restored.input_file == "/local/in"
    assert restored.output_file == "/local/stdout"
    assert restored.input_files == ("/local/in",)
    assert restored.output_files == ("/local/stdout", "/local/stderr", "/local/data")
    assert restored.manifests == ("/local/in", "/local/stdout", "/local/stderr", "/local/data")
    assert restored.all_versions == ()
    with pytest.raises(ValueError, match="availability"):
        FileEntry("/local/test", availability_state="inconsistent")


def test_execution_backend_id_is_write_once():
    execution = Execution(ExecutorEndpoint(), WorkingDirectory("/tmp"))
    submitted = execution.submitted("job-12")
    assert execution.backend_id is None
    assert submitted.backend_id == "job-12"
    assert submitted.submitted("job-12") == submitted
    with pytest.raises(ValueError, match="Cannot replace"):
        submitted.submitted("job-13")

def test_each_artifact_tracks_its_own_transfer_outcome():
    manifest = FileManifest(
        outputs=(
            FileEntry("/local/one", "/remote/one", "pending"),
            FileEntry("/local/two", "/remote/two", "pending"),
        )
    )
    manifest = manifest.with_file_availability("/local/one", "available")
    manifest = manifest.with_file_availability(
        "/local/two", "failed", "connection lost",
    )
    assert manifest.outputs[0].availability_state == "available"
    assert manifest.outputs[1].availability_state == "failed"
    assert manifest.failed_artifacts == (manifest.outputs[1],)
    assert manifest.availability_state == "failed"
    assert FileManifest(**asdict(manifest)) == manifest


def test_individual_wildcard_matches_are_recorded_separately():
    manifest = FileManifest(
        outputs=(FileEntry("/local/*.dat", "/remote/*.dat", "pending"),)
    )
    manifest = manifest.with_file_availability(
        "/local/one.dat", "available", remote_path="/remote/one.dat",
    )
    manifest = manifest.with_file_availability(
        "/local/two.dat", "failed", "I/O error", remote_path="/remote/two.dat",
    )
    assert len(manifest.outputs) == 3
    assert manifest.outputs[1].availability_state == "available"
    assert manifest.outputs[2].availability_state == "failed"
    assert FileManifest(**asdict(manifest)) == manifest

def test_partial_download_preserves_completed_files_and_cache_status():
    manifest = FileManifest(
        outputs=(
            FileEntry("/local/one", "/remote/one", "pending"),
            FileEntry("/local/two", "/remote/two", "pending"),
        ),
        caches=(CacheEntry(FileEntry("/local/cache", "/remote/cache", "pending")),),
    )
    manifest = manifest.with_file_availability("/local/one", "available")
    manifest = manifest.with_file_availability(
        "/local/two", "failed", "lost connection",
    )
    manifest = manifest.with_file_availability(
        "/local/cache", "skipped", kind="cache",
    )
    assert [x.availability_state for x in manifest.artifacts] == [
        "available", "failed", "skipped",
    ]
    assert manifest.failed_artifacts[0].local_path == "/local/two"
    assert manifest.availability_message == "lost connection"


def test_merged_stdout_stderr_has_only_one_file_status():
    manifest = FileManifest(
        stdout=FileEntry("/local/combined", "/remote/combined", "pending"),
        stderr=FileEntry("/local/combined", "/remote/combined", "pending"),
    )
    manifest = manifest.with_file_availability("/local/combined", "available")
    assert len(manifest.artifacts) == 1
    assert manifest.stdout.availability_state == "available"
    assert manifest.stderr.availability_state == "available"
    assert FileManifest(**asdict(manifest)) == manifest


def test_wildcard_parent_remains_pending_until_completion():
    manifest = FileManifest(
        outputs=(FileEntry("/local/*.dat", "/remote/*.dat", "pending"),)
    )
    manifest = manifest.with_file_availability(
        "/local/one.dat", "pending", remote_path="/remote/one.dat",
        source_pattern="/local/*.dat",
    )
    assert manifest.outputs[0].availability_state == "pending"
    assert manifest.outputs[1].source_pattern == "/local/*.dat"
    manifest = manifest.with_file_availability("/local/one.dat", "available")
    assert manifest.outputs[0].availability_state == "pending"
    manifest = manifest.with_file_availability(
        "/local/two.dat", "failed", "transfer failed",
        remote_path="/remote/two.dat", source_pattern="/local/*.dat",
    )
    manifest = manifest.with_file_availability(
        "/local/*.dat", "failed", "transfer failed",
    )
    assert [entry.availability_state for entry in manifest.outputs] == [
        "failed", "available", "failed",
    ]
    assert FileManifest(**asdict(manifest)) == manifest


def test_wildcard_parent_completed_after_all_files():
    manifest = FileManifest(
        outputs=(FileEntry("/local/*.dat", "/remote/*.dat", "pending"),)
    )
    for item in ("one", "two"):
        manifest = manifest.with_file_availability(
            f"/local/{item}.dat", "pending",
            remote_path=f"/remote/{item}.dat", source_pattern="/local/*.dat",
        )
        manifest = manifest.with_file_availability(f"/local/{item}.dat", "available")
        assert manifest.outputs[0].availability_state == "pending"
    manifest = manifest.with_file_availability("/local/*.dat", "available")
    assert manifest.availability_state == "available"
    assert len(manifest.outputs) == 3
