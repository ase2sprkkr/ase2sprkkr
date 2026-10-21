import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import (
    YamlConfigurationWarning,
    load_configuration_yaml,
    load_yaml_document,
    save_configuration_yaml,
)


def fresh_config():
    return definition.create_object()


@pytest.mark.parametrize(
    (
        "invalid_yaml", "valid_yaml", "valid_name", "expected_mpi",
        "expected_dir", "invalid_name",
    ),
    (
        (
            "  broken:\n    mpi:\n      nonsense: true\n"
            "    executable_dir: /opt/broken\n",
            "  local:\n    mpi: 4\n    executable_dir: /opt/local\n",
            "local", 4, "/opt/local", "broken",
        ),
        (
            "  broken:\n    mpi: 8\n    executble_dir: /typo\n",
            "  valid:\n    mpi: 2\n",
            "valid", 2, None, "broken",
        ),
        (
            "  ase:\n    mpi: 8\n",
            "  local:\n    mpi: 2\n",
            "local", 2, None, "ase",
        ),
    ),
    ids=("invalid-mpi", "unknown-profile-key", "reserved-ase"),
)
def test_invalid_profile_does_not_disable_other_profiles(
    tmp_path, invalid_yaml, valid_yaml, valid_name, expected_mpi,
    expected_dir, invalid_name,
):
    path = tmp_path / "config.yaml"
    path.write_text("profiles:\n" + valid_yaml + invalid_yaml, encoding="utf8")

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match=f"profile '{invalid_name}'"):
        load_configuration_yaml(config, path)

    assert config.profiles[valid_name].mpi() == expected_mpi
    if expected_dir is not None:
        assert config.profiles[valid_name].executable_dir() == expected_dir
    # ASE can still be exposed as a virtual profile from host configuration.
    if invalid_name != "ase":
        assert invalid_name not in config.profiles
    assert invalid_name not in config.profiles._values


def test_invalid_section_option_does_not_disable_siblings(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "running:\n"
        "  print_output: definitely-not-valid\n"
        "  mpi_warning: false\n"
        "nomad:\n"
        "  token: valid-token\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match="running.print_output"):
        load_configuration_yaml(config, path)

    assert config.running.mpi_warning() is False
    assert config.nomad.token() == "valid-token"
    assert config.running.print_output() == "info"


@pytest.mark.parametrize(
    ("default_name", "extra_profiles"),
    (
        ("missing", ""),
        ("broken", "  broken:\n    mpi:\n      nonsense: true\n"),
    ),
    ids=("unknown-name", "rejected-profile"),
)
def test_invalid_default_profile_is_ignored_while_valid_profiles_load(
    tmp_path, default_name, extra_profiles,
):
    path = tmp_path / "config.yaml"
    path.write_text(
        f"default_profile: {default_name}\n"
        "profiles:\n"
        + extra_profiles
        + "  local:\n    mpi: 2\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning) as caught:
        document = load_configuration_yaml(config, path)

    messages = "\n".join(str(item.message) for item in caught)
    assert "default_profile" in messages
    if default_name == "broken":
        assert "profile 'broken'" in messages
    assert config.default_profile() is None
    assert default_name not in config.profiles
    assert config.profiles["local"].mpi() == 2
    assert document["default_profile"] == default_name


@pytest.mark.parametrize("local_preexisting", (True, False), ids=("existing-local", "added-local"))
def test_save_preserves_invalid_and_unknown_raw_yaml(tmp_path, local_preexisting):
    path = tmp_path / "config.yaml"
    path.write_text(
        "# keep document comment\n"
        "profiles:\n"
        + ("  local:\n    mpi: 4\n" if local_preexisting else "")
        + "  broken:\n"
        "    # keep broken profile\n"
        "    mpi:\n"
        "      nonsense: true\n"
        "unknown_future_option:\n"
        "  value: 17\n"
        "nomad:\n"
        "  token: old\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning):
        load_configuration_yaml(config, path)

    if local_preexisting:
        assert config.profiles["local"].mpi() == 4
    else:
        config.profiles.add("local").mpi.set(4)
    config.nomad.token.set("new")
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    text = path.read_text(encoding="utf8")
    assert data["profiles"]["broken"]["mpi"] == {"nonsense": True}
    assert data["profiles"]["local"]["mpi"] == 4
    assert data["unknown_future_option"] == {"value": 17}
    assert data["nomad"]["token"] == "new"
    assert "# keep document comment" in text
    assert "# keep broken profile" in text


def test_save_preserves_only_invalid_children_of_known_section(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "running:\n"
        "  print_output: definitely-not-valid\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning):
        load_configuration_yaml(config, path)

    config.running.mpi_warning.set(False)
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    assert data["running"]["print_output"] == "definitely-not-valid"
    assert "mpi_warning" not in data["running"]
    assert data["profiles"]["legacy"]["mpi_warning"] is False


def test_invalid_default_profile_survives_unrelated_save(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "default_profile: missing  # keep invalid value\n"
        "nomad:\n"
        "  token: old\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning):
        load_configuration_yaml(config, path)

    config.nomad.token.set("new")
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    text = path.read_text(encoding="utf8")
    assert data["default_profile"] == "missing"
    assert data["nomad"]["token"] == "new"
    assert "# keep invalid value" in text


def test_explicit_valid_profile_replaces_invalid_raw_profile(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  cluster:\n"
        "    mpi:\n"
        "      nonsense: true\n"
        "    obsolete_garbage: keep-me\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning):
        load_configuration_yaml(config, path)

    profile = config.profiles.add("cluster")
    profile.mpi.set(6)
    profile.executable_dir.set("/opt/sprkkr")
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    assert data["profiles"]["cluster"] == {
        "mpi": 6,
        "executable_dir": "/opt/sprkkr",
    }
