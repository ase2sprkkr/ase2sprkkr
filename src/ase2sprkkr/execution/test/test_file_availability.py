"""Contracts shared by local project and PSI/J artifact discovery."""

from ..file_availability import discover_local_files
from ..run_model import CacheEntry, FileEntry, FileManifest


def test_discovers_outputs_and_caches_with_independent_pattern_children(tmp_path):
    (tmp_path / "out-a.dat").write_bytes(b"a")
    (tmp_path / "cache-a.pot").write_bytes(b"old")
    manifest = FileManifest(
        outputs=(FileEntry(local_path=str(tmp_path / "out-*.dat"),
                           availability_state="pending"),),
        caches=(CacheEntry(FileEntry(
            local_path=str(tmp_path / "cache-*.pot"),
            availability_state="pending",
        )),),
    )
    found = discover_local_files(manifest)
    assert len(found.outputs) == 2
    assert found.outputs[0].availability_state == "available"
    assert found.outputs[1].source_pattern == str(tmp_path / "out-*.dat")
    assert len(found.caches) == 2
    assert found.caches[0].file.availability_state == "available"
    assert found.caches[1].file.source_pattern == str(tmp_path / "cache-*.pot")
    assert discover_local_files(found) == found


def test_missing_optional_files_are_skipped_without_creating_them(tmp_path):
    missing = tmp_path / "missing.out"
    files = FileManifest(
        outputs=(FileEntry(local_path=str(missing), availability_state="pending"),),
    )
    found = discover_local_files(files)
    assert found.outputs[0].availability_state == "skipped"
    assert not missing.exists()


def test_local_supervisor_discovery_does_not_override_stream_publication(tmp_path):
    output = tmp_path / "previous.out"
    output.write_bytes(b"stale")
    files = FileManifest(
        stdout=FileEntry(local_path=str(output), availability_state="pending"),
        outputs=(FileEntry(local_path=str(tmp_path / "other.out"),
                           availability_state="pending"),),
    )
    found = discover_local_files(files, include_streams=False)
    assert found.stdout.availability_state == "pending"
    assert found.outputs[0].availability_state == "skipped"


def test_clear_outputs_preserves_inputs_and_cache(tmp_path):
    from ..file_availability import clear_declared_outputs
    input_file = tmp_path / "input.dat"
    cache_file = tmp_path / "cache.dat"
    output_file = tmp_path / "output.dat"
    other = tmp_path / "result-2.dat"
    for path in (input_file, cache_file, output_file, other):
        path.write_text("old")
    clear_declared_outputs(
        (str(tmp_path / "*.dat"),),
        protected=(input_file, cache_file),
    )
    assert input_file.read_text() == "old"
    assert cache_file.read_text() == "old"
    assert not output_file.exists()
    assert not other.exists()


def test_clear_outputs_preserves_wildcard_cache_matches(tmp_path):
    from ..file_availability import clear_declared_outputs

    cache = tmp_path / "restart-1.pot"
    output = tmp_path / "result-1.pot"
    cache.write_bytes(b"restart state")
    output.write_bytes(b"stale result")
    clear_declared_outputs(
        [str(tmp_path / "*.pot")],
        protected=[str(tmp_path / "restart-*.pot")],
    )
    assert cache.read_bytes() == b"restart state"
    assert not output.exists()
