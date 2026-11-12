"""KKRunner wip2 execution profiles keep transport, backend, and staging distinct."""

from dataclasses import replace

import pytest

from kkrunner.client import RunClient
from kkrunner.client.pool import ClientPool
from kkrunner.client.transports import LocalTransport, SshTransport
from kkrunner.execution import MpiSpec

from ase2sprkkr.sprkkr.profile import (
    SPRKKRProfile,
    client_from_profile,
    create_client_from_profile,
)


def test_local_profile_uses_direct_backend_and_reuses_client():
    profile = SPRKKRProfile(mpi="auto")
    client = client_from_profile(profile)
    assert isinstance(client, RunClient)
    assert client.transport == LocalTransport()
    assert client.backend.name == "direct"
    assert client_from_profile(replace(profile, mpi=8)) is client


def test_backend_selection_is_independent_of_ssh_transport():
    profile = SPRKKRProfile(
        executor="ssh", host="example.org",
        backend="direct", ssh_options=("-p", "2222"),
        remote_python=("python3", "-I"),
        local_root="/local", remote_root="/remote",
    )
    client = create_client_from_profile(profile)
    assert isinstance(client.transport, SshTransport)
    assert client.transport.host == "example.org"
    assert client.transport.options == ("-p", "2222")
    assert client.transport.python_command == ("python3", "-I")
    assert client.backend.name == "direct"

    managed = create_client_from_profile(replace(profile, backend="process"))
    assert managed.backend.name == "process"
    assert managed.transport == client.transport


def test_ssh_default_is_managed_and_legacy_psij_is_supported():
    ssh = SPRKKRProfile(executor="ssh", host="example.org")
    assert ssh.backend_config.name == "process"
    psij = SPRKKRProfile(executor="psij", psij_backend="slurm")
    assert psij.backend_config.name == "psij"
    assert psij.backend_config.options == (("executor_name", "slurm"), ("url", None))


def test_staging_is_per_run_not_client_identity(tmp_path):
    a = SPRKKRProfile(
        executor="ssh", host="example.org",
        local_root=tmp_path, remote_root="/scratch/a", input_sync="missing",
    )
    b = replace(a, remote_root="/scratch/b", input_sync="always")
    assert client_from_profile(a) is client_from_profile(b)
    with (tmp_path / "in.inp").open("w") as input_file:
        spec = a.create_execution("kkrscf", input_file, tmp_path / "out", tmp_path)
    assert spec.staging.local_root == str(tmp_path)
    assert spec.staging.remote_root == "/scratch/a"
    assert spec.staging.input_sync == "missing"


def test_mpi_spec_is_independent_of_backend(tmp_path):
    input_file = tmp_path / "case.inp"
    input_file.write_text("input")
    spec = SPRKKRProfile(mpi="auto").create_execution(
        "kkrscf", input_file, tmp_path / "output", tmp_path, mpi_supported=True,
    )
    assert spec.mpi == MpiSpec.from_value("auto")
    assert spec.argv == ("kkrscf",)
    assert spec.mpi_argv == ("kkrscfMPI", input_file.name)
    assert spec.mpi_uses_stdin is False

    serial = SPRKKRProfile(mpi=False).create_execution(
        "kkrscf", input_file, tmp_path / "output", tmp_path, mpi_supported=True,
    )
    assert serial.mpi == MpiSpec()
    assert serial.mpi_argv is None


def test_invalid_ssh_profile_rejected():
    with pytest.raises(ValueError, match="requires 'host'"):
        SPRKKRProfile(executor="ssh").transport
    with pytest.raises(ValueError, match="remote_command"):
        SPRKKRProfile(
            executor="ssh", host="example.org", remote_command=("custom-agent",)
        ).transport


def test_client_pool_discard_and_recreate():
    pool = ClientPool()
    profile = SPRKKRProfile()
    key = (profile.transport, profile.backend_config)
    first = pool.for_config(*key)
    assert pool.for_config(*key) is first
    assert pool.discard(key) is first
    assert pool.for_config(*key) is not first
