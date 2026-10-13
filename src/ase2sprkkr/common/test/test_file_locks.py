from ase2sprkkr.common.file_locks import FileLock


def test_nonblocking_file_lock_detects_owner(tmp_path):
    path = tmp_path / "resource.lock"
    first = FileLock(path)
    second = FileLock(path)

    assert first.acquire()
    try:
        assert not second.acquire(blocking=False)
    finally:
        first.release()

    assert second.acquire(blocking=False)
    second.release()


def test_file_lock_context_manager_releases_lock(tmp_path):
    path = tmp_path / "resource.lock"

    with FileLock(path) as lock:
        assert lock.acquired

    assert not lock.acquired
    probe = FileLock(path)
    assert probe.acquire(blocking=False)
    probe.release()
