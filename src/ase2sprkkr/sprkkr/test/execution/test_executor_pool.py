from kkrunner.backends import ExecutorPool
from kkrunner.files import FileEntry, FileManifest

from ase2sprkkr.sprkkr.profile import SPRKKRProfile, executor_from_profile


def test_unused_local_profile_fields_share_backend_and_collision_namespace(tmp_path):
    from kkrunner.backends import ExecutorEndpoint
    from kkrunner.runs.conflicts import submission_conflicts
    from kkrunner.runs.manager import RunManager
    from kkrunner.runs.records import RunStore
    pool = ExecutorPool()
    plain = SPRKKRProfile()
    unused = SPRKKRProfile(host='unused', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('unused',), psij_backend='unused', psij_url='unused')
    assert plain.endpoint == unused.endpoint == ExecutorEndpoint()
    assert pool.for_profile(plain) is pool.for_profile(unused)
    manager = RunManager(RunStore(tmp_path / 'runs'))
    first = manager.new_record(directory=tmp_path, files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "shared.out"))), endpoint=plain.endpoint, metadata={'ase2sprkkr': {'task': str('SCF').upper()}}, processing_required=True)
    second = manager.new_record(directory=tmp_path, files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "shared.out"))), endpoint=unused.endpoint, metadata={'ase2sprkkr': {'task': str('SCF').upper()}}, processing_required=True)
    assert submission_conflicts(first, second)


def test_executor_pool_reuses_backend_by_endpoint():
    pool = ExecutorPool()
    first_profile = SPRKKRProfile(executor='local', mpi=8)
    second_profile = SPRKKRProfile(executor='local', mpi=64)
    first = pool.for_profile(first_profile)
    second = pool.for_profile(second_profile)
    assert first is second


def test_ssh_profiles_with_different_staging_share_endpoint_executor():
    pool = ExecutorPool()
    first_profile = SPRKKRProfile(executor='ssh', host='cluster', local_root='/home/a', remote_root='/scratch/a', input_sync='missing')
    second_profile = SPRKKRProfile(executor='ssh', host='cluster', local_root='/home/b', remote_root='/scratch/b', input_sync='always')
    assert first_profile.endpoint == second_profile.endpoint
    assert first_profile.launch_parameters != second_profile.launch_parameters
    assert pool.for_profile(first_profile) is pool.for_profile(second_profile)


def test_executor_from_profile_uses_process_default_pool():
    profile = SPRKKRProfile(executor='local', mpi=8)
    assert executor_from_profile(profile) is executor_from_profile(profile)


def test_executor_pool_can_restart_discarded_backend():
    pool = ExecutorPool()
    profile = SPRKKRProfile(executor='local')
    first = pool.for_profile(profile)
    assert pool.discard_profile(profile) is first
    second = pool.for_profile(profile)
    assert second is not first
