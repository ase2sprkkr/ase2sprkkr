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


def test_project_manifest_reuses_optional_file_path_semantics(tmp_path):
    from types import SimpleNamespace

    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file="input.dat",
        input_files=(object(), 1, "input.dat"),
        output_file="stdout.dat",
        error_file=None,
        output_files=("result.dat",),
        cache_files=("cache.dat",),
    )
    manifest = FileManifest.from_spec(spec)

    assert manifest.input_files == (str(tmp_path / "input.dat"),)
    assert manifest.stdout.local_path == str(tmp_path / "stdout.dat")
    assert manifest.stderr is None
    assert [item.local_path for item in manifest.outputs] == [
        str(tmp_path / "result.dat")
    ]
    assert manifest.caches[0].file.local_path == str(tmp_path / "cache.dat")



def test_recovered_stream_paths_update_only_declared_stream_locations():
    manifest = FileManifest(
        inputs=(FileEntry("/local/input", "/remote/input"),),
        outputs=(FileEntry("/local/artifact", "/remote/artifact"),),
        caches=(CacheEntry(FileEntry("/local/cache", "/remote/cache"), (12, 34)),),
        stdout=FileEntry(
            "/local/stdout", "/previous/stdout",
            availability_state="pending",
        ),
        stderr=FileEntry("/local/stderr"),
    )

    updated = manifest.with_remote_stream_paths(
        stdout="/recovered/out", stderr="/recovered/err",
    )

    assert updated.stdout.remote_path == "/recovered/out"
    assert updated.stderr.remote_path == "/recovered/err"
    assert updated.stdout.local_path == manifest.stdout.local_path
    assert updated.stdout.availability_state == "pending"
    assert updated.inputs == manifest.inputs
    assert updated.outputs == manifest.outputs
    assert updated.caches == manifest.caches
    assert manifest.stdout.remote_path == "/previous/stdout"
    assert updated.with_remote_stream_paths(stdout="", stderr=None) is updated


def test_recovered_stream_paths_do_not_create_unrecorded_stream():
    manifest = FileManifest(stdout=FileEntry("/local/stdout"))

    updated = manifest.with_remote_stream_paths(
        stdout="/remote/out", stderr="/remote/unrecorded.err",
    )
    assert updated.stdout.remote_path == "/remote/out"
    assert updated.stderr is None


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize(
    "state, message, error",
    [
        ("invalid", None, "availability state"),
        ("available", "unexpected error", "availability error"),
    ],
)
def test_manifest_availability_validation_is_shared_with_file_entry(existing, state, message, error):
    manifest = FileManifest(
        outputs=(FileEntry("/local/result"),) if existing else (),
    )
    with pytest.raises(ValueError, match=error):
        manifest.with_file_availability("/local/result", state, message)


def test_submitted_uses_execution_backend_id_validation():
    execution = Execution(ExecutorEndpoint(), WorkingDirectory("/tmp"))
    with pytest.raises(ValueError, match="backend_id cannot be empty"):
        execution.submitted("")
    bound = execution.submitted(42)
    assert bound.backend_id == "42"
    assert bound.submitted(42) == bound
    with pytest.raises(ValueError, match="backend_id cannot be empty"):
        bound.submitted("")
    with pytest.raises(ValueError, match="Cannot replace"):
        bound.submitted("other")


def test_manifest_from_spec_normalizes_stream_overlaps(tmp_path):
    from types import SimpleNamespace

    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file="input.inp", input_files=("input.inp", "other.inp", "other.inp"),
        output_file="stdout.out", error_file="stderr.err",
        output_files=("stdout.out", "stderr.err", "result.dat", "result.dat"),
        cache_files=("restart.pot", "restart.pot"),
    )
    files = FileManifest.from_spec(spec)
    assert files.input_files == (str(tmp_path / "input.inp"), str(tmp_path / "other.inp"))
    assert files.stdout.local_path == str(tmp_path / "stdout.out")
    assert files.stderr.local_path == str(tmp_path / "stderr.err")
    assert tuple(x.local_path for x in files.outputs) == (str(tmp_path / "result.dat"),)
    assert files.output_files == (
        str(tmp_path / "stdout.out"),
        str(tmp_path / "stderr.err"),
        str(tmp_path / "result.dat"),
    )
    assert len(files.caches) == 1
    assert FileManifest(**asdict(files)) == files


def test_legacy_duplicate_stderr_declaration_is_updated_together():
    path = "/local/stderr"
    files = FileManifest(
        stderr=FileEntry(path, availability_state="pending"),
        outputs=(FileEntry(path, availability_state="pending"),),
    )
    updated = files.with_file_availability(path, "available")
    assert updated.stderr.availability_state == "available"
    assert updated.outputs[0].availability_state == "available"
    assert updated.output_files == (path,)
