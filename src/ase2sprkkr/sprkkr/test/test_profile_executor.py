"""Named profiles compose KKRunner transport, backend, and per-run staging."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr.profile import (
    SPRKKRProfile, executor_from_profile, resolve_sprkkr_profile,
)
from kkrunner.client import EndpointClient
from kkrunner.client.transports import LocalTransport, SshTransport


def fresh_config(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config", staticmethod(lambda: None),
    )
    return definition.create_object()


def test_named_ssh_profile_resolves_client_and_staging(monkeypatch):
    configuration = fresh_config(monkeypatch)
    configuration.profiles.set({
        "cluster": {
            "executor": "ssh", "host": "user@cluster",
            "backend": "psij", "psij_backend": "slurm",
            "local_root": "/home/user/work",
            "remote_root": "/scratch/user/work",
            "ssh_options": ["-p", "2222"], "ssh_zstd": True,
            "remote_python": ["python3"],
            "input_sync": "missing", "output_sync": "always", "mpi": 16,
        },
    })
    profile = resolve_sprkkr_profile("cluster", configuration=configuration)
    assert profile.executor == "ssh"
    assert profile.host == "user@cluster"
    assert profile.backend == "psij"
    assert profile.local_root == "/home/user/work"
    assert profile.remote_root == "/scratch/user/work"
    assert profile.ssh_options == ("-p", "2222")
    assert profile.ssh_zstd is True
    assert profile.remote_python == ("python3",)
    assert profile.input_sync == "missing"
    assert profile.output_sync == "always"
    assert profile.mpi == 16

    client = executor_from_profile(profile)
    assert isinstance(client, EndpointClient)
    assert isinstance(client.transport, SshTransport)
    assert client.backend.name == "psij"


def test_local_profile_uses_direct_client():
    client = executor_from_profile(SPRKKRProfile())
    assert isinstance(client, EndpointClient)
    assert isinstance(client.transport, LocalTransport)
    assert client.backend.name == "direct"


def test_ssh_profile_does_not_embed_staging_into_transport(tmp_path):
    profile = SPRKKRProfile(
        executor="ssh", host="user@cluster",
        local_root="/home/user/work", remote_root="/scratch/user/work",
        ssh_options=("-p", "2222"), ssh_zstd=True,
        input_sync="missing", output_sync="always",
    )
    client = executor_from_profile(profile)
    assert client.transport == SshTransport(
        "user@cluster", options=("-p", "2222"), zstd=True
    )
    input_file = SimpleNamespace(name="/home/user/work/test.inp")
    spec = profile.create_execution(
        "kkrscf", input_file, "/home/user/work/out",
        "/home/user/work", mpi_supported=True,
    )
    assert spec.staging.local_root == "/home/user/work"
    assert spec.staging.remote_root == "/scratch/user/work"
    assert spec.staging.input_sync == "missing"
    assert spec.staging.output_sync == "always"


def test_ssh_profile_requires_host():
    with pytest.raises(ValueError, match="requires 'host'"):
        executor_from_profile(SPRKKRProfile(executor="ssh"))


def test_psij_profile_creates_configured_backend():
    client = executor_from_profile(SPRKKRProfile(
        executor="psij", psij_backend="slurm", psij_url="local://scheduler",
    ))
    assert isinstance(client.transport, LocalTransport)
    assert client.backend.name == "psij"
    assert client.backend.options == {
        "executor_name": "slurm", "url": "local://scheduler",
    }


def test_unknown_runtime_executor_is_rejected():
    with pytest.raises(ValueError, match="Unknown execution location"):
        executor_from_profile(SPRKKRProfile(executor="custom"))


def test_mpi_intent_does_not_change_pooled_client(tmp_path):
    profile = SPRKKRProfile(mpi=8)
    one = executor_from_profile(profile)
    other = executor_from_profile(replace(profile, mpi=16))
    assert one is other
    path = tmp_path / "run.inp"
    path.write_text("input")
    spec = profile.create_execution(
        "kkrscf", path, tmp_path / "run.out", tmp_path, mpi_supported=True,
    )
    assert spec.mpi.required
    assert spec.mpi.process_count == 8
