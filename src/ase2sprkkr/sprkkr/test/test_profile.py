import shlex

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import load_configuration_yaml
from ase2sprkkr.sprkkr.profile import (
    SPRKKRProfile,
    profile_from_ase_command,
    resolve_sprkkr_profile,
)


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
