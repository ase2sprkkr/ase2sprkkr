from kkrunner.backends.config import BackendConfig
from kkrunner.client.transports import SshTransport

from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def test_profile_persists_transport_separately_from_backend():
    profile = SPRKKRProfile(
        executor="ssh", host="cluster",
        ssh_options=("-p", "2222"), ssh_zstd=True,
        remote_python=("python3",),
        backend="psij", psij_backend="slurm", psij_url="local://slurm",
    )
    assert profile.transport == SshTransport(
        host="cluster", options=("-p", "2222"), zstd=True,
        python_command=("python3",),
    )
    assert profile.backend_config == BackendConfig(
        "psij", {"executor_name": "slurm", "url": "local://slurm"},
    )
