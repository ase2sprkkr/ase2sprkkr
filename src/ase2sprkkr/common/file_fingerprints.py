"""Cheap filesystem fingerprints for detecting replaced or rewritten files."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple


class FileFingerprint(NamedTuple):
    """Metadata identity of one filesystem object at one point in time.

    This intentionally avoids reading file contents.  It is suitable for detecting
    ordinary replacement or rewrite races around potentially large SPR-KKR files,
    not for cryptographic content identity.
    """

    device: int
    inode: int
    mode: int
    size: int
    mtime_ns: int
    ctime_ns: int


def file_fingerprint(path: str | Path) -> FileFingerprint:
    """Return a cheap fingerprint for ``path`` based on one ``stat()`` call."""
    info = Path(path).stat()
    return FileFingerprint(
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )
