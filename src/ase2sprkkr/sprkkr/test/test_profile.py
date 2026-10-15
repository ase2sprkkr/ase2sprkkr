import shlex

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import load_configuration_yaml
from ase2sprkkr.sprkkr.profile import (
    ExecutorEndpoint,
    LaunchParameters,
    SPRKKRProfile,
    profile_from_ase_command,
    resolve_sprkkr_profile,
)


def test_configured_profiles_expose_mutability_contract():
    config = definition.create_object()

    profile = config.profiles.add("cluster")

    assert profile.read_only is False
    assert config.profiles["cluster"].read_only is False


def test_profile_shell_strings_round_trip_to_runtime_argv(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  cluster:\n"
        "    mpi: srun -n 8\n"
        "    ssh_options: -p 2222 -i '/tmp/key file'\n",
        encoding="utf8",
    )
    config = definition.create_object()

    load_configuration_yaml(config, path)
    profile = resolve_sprkkr_profile("cluster", configuration=config)

    assert profile.mpi == ("srun", "-n", "8")
    assert profile.ssh_options == ("-p", "2222", "-i", "/tmp/key file")


def test_explicit_runtime_shell_strings_are_split():
    profile = resolve_sprkkr_profile(
        SPRKKRProfile(
            mpi="mpirun -np 4 --bind-to 'core list'",
            ssh_options="-o 'ProxyCommand=ssh jump nc %h %p' -p 2222",
        )
    )

    assert profile.mpi == tuple(shlex.split("mpirun -np 4 --bind-to 'core list'"))
    assert profile.ssh_options == tuple(
        shlex.split("-o 'ProxyCommand=ssh jump nc %h %p' -p 2222")
    )


def test_mpi_auto_keeps_symbolic_semantics():
    profile = resolve_sprkkr_profile(SPRKKRProfile(mpi="auto"))

    assert profile.mpi == "auto"


def test_ase_command_locates_sprkkr_executable_family():
    profile = profile_from_ase_command("/opt/sprkkr/kkrscf_v8")

    assert profile.executable_dir == "/opt/sprkkr"
    assert profile.executable_suffix == "_v8"


def test_runtime_profile_splits_endpoint_from_launch_parameters():
    profile = SPRKKRProfile(
        executor="ssh",
        host="cluster.example",
        ssh_options=("-p", "2222"),
        remote_command=("python", "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
        psij_backend="slurm",
        psij_url="ssh://scheduler",
        executable_dir="/opt/sprkkr",
        executable_suffix="_v8",
        mpi=32,
        local_root="/home/me/project",
        remote_root="/scratch/me/project",
        input_sync="missing",
        output_sync="always",
    )

    assert profile.endpoint == ExecutorEndpoint(
        executor="ssh",
        host="cluster.example",
        ssh_options=("-p", "2222"),
        remote_command=("python", "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
        psij_backend="slurm",
        psij_url="ssh://scheduler",
    )
    assert profile.launch_parameters == LaunchParameters(
        executable_dir="/opt/sprkkr",
        executable_suffix="_v8",
        mpi=32,
        local_root="/home/me/project",
        remote_root="/scratch/me/project",
        input_sync="missing",
        output_sync="always",
    )


def test_profiles_with_same_backend_share_endpoint_identity():
    first = SPRKKRProfile(
        executor="ssh",
        host="cluster.example",
        psij_backend="slurm",
        mpi=8,
        remote_root="/scratch/one",
    )
    second = SPRKKRProfile(
        executor="ssh",
        host="cluster.example",
        psij_backend="slurm",
        mpi=64,
        remote_root="/scratch/two",
        input_sync="never",
    )

    assert first.endpoint == second.endpoint
    assert first.launch_parameters != second.launch_parameters


@pytest.mark.parametrize(
    "command",
    (
        "mpirun -n 8 /opt/sprkkr/kkrscf",
        "srun /opt/sprkkr/kkrscf",
        "/opt/sprkkr/kkrspec",
        "/opt/sprkkr/{}.mpi",
    ),
)
def test_ase_command_requires_direct_kkrscf(command):
    with pytest.raises(ValueError, match="kkrscf|prefixes|wrappers"):
        profile_from_ase_command(command)
