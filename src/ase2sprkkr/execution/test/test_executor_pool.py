import pytest

from ase2sprkkr.execution import ExecutorPool
from ase2sprkkr.sprkkr.profile import SPRKKRProfile, executor_from_profile


def test_executor_pool_reuses_backend_by_endpoint():
    pool = ExecutorPool()
    first_profile = SPRKKRProfile(executor="local", mpi=8)
    second_profile = SPRKKRProfile(executor="local", mpi=64)

    first = pool.for_profile(first_profile)
    second = pool.for_profile(second_profile)

    assert first is second


def test_ssh_profiles_with_different_staging_share_endpoint_executor():
    pool = ExecutorPool()
    first_profile = SPRKKRProfile(
        executor="ssh",
        host="cluster",
        local_root="/home/a",
        remote_root="/scratch/a",
        input_sync="missing",
    )
    second_profile = SPRKKRProfile(
        executor="ssh",
        host="cluster",
        local_root="/home/b",
        remote_root="/scratch/b",
        input_sync="always",
    )

    assert first_profile.endpoint == second_profile.endpoint
    assert first_profile.launch_parameters != second_profile.launch_parameters
    assert pool.for_profile(first_profile) is pool.for_profile(second_profile)


def test_executor_from_profile_uses_process_default_pool():
    profile = SPRKKRProfile(executor="local", mpi=8)

    assert executor_from_profile(profile) is executor_from_profile(profile)


def test_executor_pool_can_restart_discarded_backend():
    pool = ExecutorPool()
    profile = SPRKKRProfile(executor="local")

    first = pool.for_profile(profile)
    assert pool.discard_profile(profile) is first
    second = pool.for_profile(profile)

    assert second is not first


def test_executor_pool_close_attempts_all_after_failure():
    """One failing backend must not prevent cleanup of the remaining backends."""
    pool = ExecutorPool()
    closed = []
    backends = {}

    class Executor:
        def __init__(self, key, error=None):
            self.key = key
            self.error = error

        def close(self):
            closed.append(self.key)
            if self.error is not None:
                raise self.error

    first_error = OSError("first close failed")
    for key, error in (
        ("first", first_error),
        ("second", RuntimeError("second close failed")),
        ("third", None),
    ):
        backend = Executor(key, error)
        backends[key] = backend
        assert pool.get(key, lambda backend=backend: backend) is backend

    with pytest.raises(OSError, match="first close failed") as caught:
        pool.close()

    assert caught.value is first_error
    assert closed == ["first", "second", "third"]
    replacement = object()
    assert pool.get("first", lambda: replacement) is replacement
    assert replacement is not backends["first"]
    pool.close()
    assert closed == ["first", "second", "third"]


def test_executor_pool_close_accepts_backends_without_close():
    pool = ExecutorPool()
    executor = object()
    pool.get("no-close", lambda: executor)

    pool.close()

    assert pool.get("no-close", object) is not executor
