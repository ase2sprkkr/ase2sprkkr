from __future__ import annotations

import os

from ase2sprkkr.common.file_versions import FileState, FileVersion


def test_metadata_only_change_keeps_version_current(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"same contents")
    version = FileVersion.capture(path)
    original = path.stat()

    os.utime(
        path,
        ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000_000),
    )

    assert version.digest is not None
    assert version.current() == version
    assert version.state() is FileState.CURRENT


def test_identical_replacement_keeps_version_current(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"same contents")
    version = FileVersion.capture(path)

    replacement = tmp_path / "replacement.dat"
    replacement.write_bytes(b"same contents")
    os.replace(replacement, path)

    assert version.current() == version
    assert version.state() is FileState.CURRENT


def test_same_size_rewrite_with_restored_mtime_is_overwritten(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"first")
    version = FileVersion.capture(path)
    original = path.stat()

    path.write_bytes(b"other")
    os.utime(
        path,
        ns=(original.st_atime_ns, original.st_mtime_ns),
    )

    assert version.state() is FileState.OVERWRITTEN
    assert version.current() != version


def test_missing_version_reports_missing(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"contents")
    version = FileVersion.capture(path)
    path.unlink()

    assert version.state() is FileState.MISSING


def test_metadata_only_snapshot_is_not_available_content_identity(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"contents")

    version = FileVersion.capture(path, digest=False)

    assert version.digest is None
    assert version.available is False
