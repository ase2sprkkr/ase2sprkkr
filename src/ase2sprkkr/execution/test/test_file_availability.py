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


def test_local_wildcard_discovery_uses_same_glob_semantics_as_transport(tmp_path):
    from ..transport import expand_local_candidates
    hidden = tmp_path / ".hidden"
    hidden.mkdir()
    (hidden / "invisible.dat").write_text("hidden")
    shown = tmp_path / "shown"
    shown.mkdir()
    (shown / "file.dat").write_text("visible")

    declaration = str(tmp_path / "**" / "*.dat")
    assert expand_local_candidates(declaration) == (str(shown / "file.dat"),)
    manifest = FileManifest(
        outputs=(FileEntry(declaration, availability_state="pending"),),
    )
    discovered = discover_local_files(manifest)
    assert [entry.local_path for entry in discovered.outputs] == [
        declaration, str(shown / "file.dat")
    ]


def test_cleanup_unlinks_symlink_without_deleting_its_target(tmp_path):
    from ..file_availability import clear_declared_outputs
    outside = tmp_path / "outside"
    outside.mkdir()
    original = outside / "original.dat"
    original.write_text("important")
    working = tmp_path / "working"
    working.mkdir()
    link = working / "result.dat"
    link.symlink_to(original)
    clear_declared_outputs((str(working / "*.dat"),))
    assert not link.is_symlink()
    assert original.read_text() == "important"


def test_protected_output_stream_matches_cache_glob(tmp_path):
    import pytest
    from ase2sprkkr.execution.file_availability import validate_output_collisions

    protected = tmp_path / "restart-1.pot"
    protected.write_text("density")
    with pytest.raises(ValueError, match="overlaps protected"):
        validate_output_collisions(
            (protected,), protected=(tmp_path / "restart-*.pot",),
        )


def test_explicit_output_collision_preserves_protected_file(tmp_path):
    import pytest
    from ase2sprkkr.execution.file_availability import validate_output_collisions

    input_path = tmp_path / "input.dat"
    input_path.write_bytes(b"input")
    with pytest.raises(ValueError, match="overlaps protected"):
        validate_output_collisions((input_path,), protected=(input_path,))
    assert input_path.read_bytes() == b"input"


def test_output_pattern_can_skip_protected_files(tmp_path):
    from ase2sprkkr.execution.file_availability import (
        validate_output_collisions, clear_declared_outputs,
    )
    keep = tmp_path / "restart.pot"
    drop = tmp_path / "old.pot"
    keep.write_text("restart")
    drop.write_text("old")
    validate_output_collisions((tmp_path / "*.pot",), protected=(keep,))
    clear_declared_outputs((tmp_path / "*.pot",), protected=(keep,))
    assert keep.read_text() == "restart"
    assert not drop.exists()


def test_protected_output_rejects_symlink_alias(tmp_path):
    import pytest
    from ase2sprkkr.execution.file_availability import validate_output_collisions

    source = tmp_path / "source.dat"
    source.write_bytes(b"density")
    alias = tmp_path / "alias.out"
    alias.symlink_to(source)
    with pytest.raises(ValueError, match="overlaps protected"):
        validate_output_collisions((alias,), protected=(source,))
    assert source.read_bytes() == b"density"


def test_protected_output_rejects_hardlink_alias(tmp_path):
    import os
    import pytest
    from ase2sprkkr.execution.file_availability import validate_output_collisions

    source = tmp_path / "source.dat"
    source.write_bytes(b"density")
    alias = tmp_path / "alias.out"
    os.link(source, alias)
    with pytest.raises(ValueError, match="overlaps protected"):
        validate_output_collisions((alias,), protected=(source,))
    assert source.read_bytes() == b"density"
