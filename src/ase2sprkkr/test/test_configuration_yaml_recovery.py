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


def test_invalid_profile_does_not_disable_other_profiles(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  local:\n"
        "    mpi: 4\n"
        "    executable_dir: /opt/local\n"
        "  broken:\n"
        "    mpi:\n"
        "      nonsense: true\n"
        "    executable_dir: /opt/broken\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match="profile 'broken'"):
        load_configuration_yaml(config, path)

    assert config.profiles["local"].mpi() == 4
    assert config.profiles["local"].executable_dir() == "/opt/local"
    assert "broken" not in config.profiles


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


def test_unknown_profile_member_invalidates_whole_profile(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  broken:\n"
        "    mpi: 8\n"
        "    executble_dir: /typo\n"
        "  valid:\n"
        "    mpi: 2\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match="profile 'broken'"):
        load_configuration_yaml(config, path)

    assert "broken" not in config.profiles
    assert config.profiles["valid"].mpi() == 2


def test_stored_ase_profile_is_rejected_but_other_profiles_load(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  ase:\n"
        "    mpi: 8\n"
        "  local:\n"
        "    mpi: 2\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match="profile 'ase'"):
        load_configuration_yaml(config, path)

    assert config.profiles["local"].mpi() == 2
    assert "ase" not in config.profiles._values


def test_unknown_default_profile_is_ignored(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "default_profile: missing\n"
        "profiles:\n"
        "  local:\n"
        "    mpi: 2\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning, match="default_profile"):
        document = load_configuration_yaml(config, path)

    assert config.default_profile() is None
    assert config.profiles["local"].mpi() == 2
    assert document["default_profile"] == "missing"


def test_default_to_rejected_profile_is_ignored(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "default_profile: broken\n"
        "profiles:\n"
        "  broken:\n"
        "    mpi:\n"
        "      nonsense: true\n"
        "  local:\n"
        "    mpi: 2\n",
        encoding="utf8",
    )

    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning) as caught:
        load_configuration_yaml(config, path)

    message = str(caught[0].message)
    assert "profile 'broken'" in message
    assert "default_profile" in message
    assert config.default_profile() is None
    assert "broken" not in config.profiles
    assert config.profiles["local"].mpi() == 2


def test_save_preserves_invalid_and_unknown_raw_yaml(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "# keep document comment\n"
        "profiles:\n"
        "  local:\n"
        "    mpi: 4\n"
        "  broken:\n"
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

    config.nomad.token.set("new")
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    text = path.read_text(encoding="utf8")
    assert data["profiles"]["broken"]["mpi"] == {"nonsense": True}
    assert data["unknown_future_option"] == {"value": 17}
    assert data["nomad"]["token"] == "new"
    assert "# keep document comment" in text
    assert "# keep broken profile" in text


def test_save_preserves_only_invalid_profile_when_adding_valid_profile(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "profiles:\n"
        "  broken:\n"
        "    mpi:\n"
        "      nonsense: true\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    with pytest.warns(YamlConfigurationWarning):
        load_configuration_yaml(config, path)

    config.profiles.add("local").mpi.set(4)
    save_configuration_yaml(config, path, document)

    data = load_yaml_document(path)
    assert data["profiles"]["broken"]["mpi"] == {"nonsense": True}
    assert data["profiles"]["local"]["mpi"] == 4


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
    assert data["running"]["mpi_warning"] is False


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
