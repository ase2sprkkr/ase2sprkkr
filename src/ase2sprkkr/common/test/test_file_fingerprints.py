from __future__ import annotations

from ase2sprkkr.common.file_fingerprints import FileFingerprint, file_fingerprint


def test_file_fingerprint_matches_stat_and_changes_on_rewrite(tmp_path):
    path = tmp_path / "artifact.dat"
    path.write_bytes(b"old")

    fingerprint = file_fingerprint(path)
    info = path.stat()

    assert fingerprint == FileFingerprint(
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )

    path.write_bytes(b"replacement is longer")
    assert file_fingerprint(path) != fingerprint
