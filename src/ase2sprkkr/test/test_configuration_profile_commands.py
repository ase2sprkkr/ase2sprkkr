from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import load_yaml_document, save_configuration_yaml


def test_profile_commands_are_saved_as_shell_strings(tmp_path):
    path = tmp_path / "config.yaml"
    config = definition.create_object()
    profile = config.profiles.add("cluster")
    profile.mpi.set(["srun", "-n", "8"])
    profile.ssh_options.set(["-p", "2222", "-i", "/tmp/key file"])

    save_configuration_yaml(config, path)
    data = load_yaml_document(path)

    assert data["profiles"]["cluster"]["mpi"] == "srun -n 8"
    assert data["profiles"]["cluster"]["ssh_options"] == (
        "-p 2222 -i '/tmp/key file'"
    )
