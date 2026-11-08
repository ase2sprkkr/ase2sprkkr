from kkrunner.backends import ExecutorEndpoint

from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def test_profile_uses_execution_endpoint_value():
    profile = SPRKKRProfile(executor='ssh', host='cluster', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('kkrunner', 'remote'), psij_backend='slurm', psij_url='local://slurm')
    assert profile.endpoint == ExecutorEndpoint(executor='ssh', host='cluster', ssh_options=('-p', '2222'), ssh_zstd=True, remote_command=('kkrunner', 'remote'), psij_backend='slurm', psij_url='local://slurm')
