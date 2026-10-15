"""Cross-platform inter-process file locking helpers."""

from __future__ import annotations

import os
from pathlib import Path


class FileLock:
    """Exclusive inter-process lock backed by a one-byte lock file."""

    def __init__(self, path):
        self.path = Path(path)
        self._file = None

    @property
    def acquired(self):
        return self._file is not None

    @staticmethod
    def _ensure_lock_byte(file):
        file.seek(0, os.SEEK_END)
        if file.tell() == 0:
            file.write(b"\0")
            file.flush()
        file.seek(0)

    def acquire(self, *, blocking=True):
        """Acquire the lock and return whether it was acquired."""
        if self._file is not None:
            raise RuntimeError(f"File lock is already acquired: {self.path}")

        self.path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        file = self.path.open("a+b")
        try:
            if os.name == "nt":
                import msvcrt

                self._ensure_lock_byte(file)
                mode = msvcrt.LK_LOCK if blocking else msvcrt.LK_NBLCK
                try:
                    msvcrt.locking(file.fileno(), mode, 1)
                except OSError:
                    if not blocking:
                        file.close()
                        return False
                    raise
            else:
                import fcntl

                operation = fcntl.LOCK_EX
                if not blocking:
                    operation |= fcntl.LOCK_NB
                try:
                    fcntl.flock(file.fileno(), operation)
                except BlockingIOError:
                    file.close()
                    return False
        except BaseException:
            if not file.closed:
                file.close()
            raise

        self._file = file
        return True

    def release(self):
        """Release the lock. Releasing an unlocked instance is harmless."""
        file = self._file
        if file is None:
            return
        self._file = None
        try:
            if os.name == "nt":
                import msvcrt

                file.seek(0)
                msvcrt.locking(file.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(file.fileno(), fcntl.LOCK_UN)
        finally:
            file.close()

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.release()
