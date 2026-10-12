import shlex

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import load_configuration_yaml
from ase2sprkkr.sprkkr.profile import SPRKKRProfile, resolve_sprkkr_profile


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
