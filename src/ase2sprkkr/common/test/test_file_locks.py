from ase2sprkkr.common.file_locks import FileLock


def test_nonblocking_file_lock_ownership_and_context_release(tmp_path):
    path = tmp_path / "resource.lock"
    first = FileLock(path)
    contender = FileLock(path)

    # A second owner must not acquire the locked file.
    with first as acquired:
        assert acquired.acquired
        assert not contender.acquire(blocking=False)

    # Context exit releases the lock; the contender can then acquire it.
    assert not first.acquired
    assert contender.acquire(blocking=False)
    contender.release()

    # Explicit release must also make the same resource reacquirable.
    assert first.acquire(blocking=False)
    try:
        assert not contender.acquire(blocking=False)
    finally:
        first.release()
    assert contender.acquire(blocking=False)
    contender.release()
