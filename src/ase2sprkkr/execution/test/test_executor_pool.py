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
