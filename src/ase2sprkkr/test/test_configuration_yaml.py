import pytest

import ase2sprkkr.configuration as configuration
from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_migration import (
    clean_legacy_configuration,
    migrate_legacy_configuration,
)
from ase2sprkkr.configuration_yaml import (
    YamlConfigurationWarning,
    load_configuration_yaml,
    load_yaml_document,
    save_configuration_yaml,
)


def fresh_config():
    return definition.create_object()


def test_yaml_round_trip_named_profiles(tmp_path):
    path = tmp_path / "config.yaml"
    config = fresh_config()
    profile = config.profiles.add("local")
    profile.mpi.set(4)
    profile.mpi_warning.set(False)
    profile.executable_dir.set("/opt/sprkkr")
    config.default_profile.set("local")

    save_configuration_yaml(config, path)

    loaded = fresh_config()
    load_configuration_yaml(loaded, path)

    assert loaded.default_profile() == "local"
    assert loaded.profiles["local"].mpi() == 4
    assert loaded.profiles["local"].mpi_warning() is False
    assert loaded.profiles["local"].executable_dir() == "/opt/sprkkr"


def test_yaml_does_not_persist_implicit_profile_mpi_defaults(tmp_path):
    path = tmp_path / "config.yaml"
    config = fresh_config()
    profile = config.profiles.add("local")

    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is True
    save_configuration_yaml(config, path)

    data = load_yaml_document(path)
    stored = (data.get("profiles") or {}).get("local", {})
    assert "mpi" not in stored
    assert "mpi_warning" not in stored


def test_yaml_persists_legacy_backing_profile(tmp_path):
    path = tmp_path / "config.yaml"
    config = fresh_config()
    config.running.mpi.set(8)
    config.running.mpi_warning.set(False)
    config.executables.dir.set("/opt/legacy")

    save_configuration_yaml(config, path)
    data = load_yaml_document(path)

    assert data["profiles"]["legacy"] == {
        "mpi": 8,
        "mpi_warning": False,
        "executable_dir": "/opt/legacy",
    }
    assert "running" not in data or "mpi" not in data["running"]
    assert "running" not in data or "mpi_warning" not in data["running"]
    assert "executables" not in data or "dir" not in data["executables"]

    loaded = fresh_config()
    load_configuration_yaml(loaded, path)
    assert loaded.running.mpi() == 8
    assert loaded.running.mpi_warning() is False
    assert loaded.executables.dir() == "/opt/legacy"
    with pytest.raises(TypeError, match="read-only"):
        loaded.profiles["legacy"].mpi.set(2)


def test_set_permanent_writes_typed_yaml(tmp_path, monkeypatch):
    path = tmp_path / "config.yaml"
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(path))
    config = fresh_config()

    config.running.mpi.set_permanent(6)
    config.nomad.token.set_permanent("secret")

    data = load_yaml_document(path)
    assert data["profiles"]["legacy"]["mpi"] == 6
    assert isinstance(data["profiles"]["legacy"]["mpi"], int)
    assert data["nomad"]["token"] == "secret"


def test_set_permanent_preserves_existing_yaml_comments(tmp_path, monkeypatch):
    path = tmp_path / "config.yaml"
    path.write_text(
        "# keep me\n"
        "nomad:\n"
        "  token: old  # token comment\n",
        encoding="utf8",
    )
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(path))
    config = fresh_config()
    load_configuration_yaml(config, path)

    config.nomad.token.set_permanent("new")

    text = path.read_text(encoding="utf8")
    assert "# keep me" in text
    assert "# token comment" in text
    assert "token: new" in text


def test_yaml_does_not_persist_ase_profile(tmp_path, monkeypatch):
    import configparser
    from types import SimpleNamespace
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"command": "sprkkr"}
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )

    path = tmp_path / "config.yaml"
    config = fresh_config()
    assert "ase" in config.profiles

    save_configuration_yaml(config, path)
    data = load_yaml_document(path)

    assert "ase" not in (data.get("profiles") or {})


def test_yaml_save_preserves_comments(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "# User configuration\n"
        "default_profile: local  # preferred machine\n"
        "profiles:\n"
        "  local:\n"
        "    # executable installation\n"
        "    executable_dir: /old/path\n",
        encoding="utf8",
    )

    document = load_yaml_document(path)
    config = fresh_config()
    load_configuration_yaml(config, path)
    config.profiles["local"].executable_dir.set("/new/path")
    save_configuration_yaml(config, path, document)

    text = path.read_text(encoding="utf8")
    assert "# User configuration" in text
    assert "# preferred machine" in text
    assert "# executable installation" in text
    assert "executable_dir: /new/path" in text


def test_yaml_requires_mapping_root(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("- one\n- two\n", encoding="utf8")

    with pytest.raises(ValueError, match="must contain a mapping"):
        load_yaml_document(path)


def test_yaml_reports_unknown_configuration_keys(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "unknown_option: value\n"
        "nomad:\n"
        "  token: still-loaded\n",
        encoding="utf8",
    )
    config = fresh_config()

    with pytest.warns(YamlConfigurationWarning, match="unknown_option"):
        load_configuration_yaml(config, path)

    assert config.nomad.token() == "still-loaded"


def test_user_preferences_prefer_yaml_over_legacy_python(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"
    yaml_path.write_text("default_profile: yaml\nprofiles:\n  yaml: {}\n", encoding="utf8")
    legacy_path.write_text("config.default_profile = 'python'\n", encoding="utf8")

    config = fresh_config()
    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(configuration, "legacy_user_preferences_file", lambda: str(legacy_path))

    configuration.load_user_preferences()

    assert config.default_profile() == "yaml"
    assert "yaml" in config.profiles


def test_user_preferences_fall_back_to_legacy_python(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"
    legacy_path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.default_profile = 'legacy'\n"
        "config.profiles.add('legacy-name')\n",
        encoding="utf8",
    )

    config = fresh_config()
    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(configuration, "legacy_user_preferences_file", lambda: str(legacy_path))

    configuration.load_user_preferences()

    assert config.default_profile() == "legacy"
    assert "legacy-name" in config.profiles


def test_legacy_python_generated_options_load_into_legacy_profile(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"
    legacy_path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.executables.suffix = '_legacy'\n"
        "config.executables.dir = '/opt/sprkkr'\n"
        "config.running.mpi = 4\n",
        encoding="utf8",
    )

    config = fresh_config()
    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(configuration, "legacy_user_preferences_file", lambda: str(legacy_path))

    configuration.load_user_preferences()

    assert config.executables.suffix() == "_legacy"
    assert config.executables.dir() == "/opt/sprkkr"
    assert config.running.mpi() == 4
    assert config.running.mpi_warning() is True
    assert config.profiles["legacy"].as_dict() == {
        "mpi": 4,
        "mpi_warning": True,
        "executable_dir": "/opt/sprkkr",
        "executable_suffix": "_legacy",
    }


def test_legacy_migration_is_self_contained(tmp_path):
    legacy_path = tmp_path / "__init__.py"
    yaml_path = tmp_path / "config.yaml"
    legacy_path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.executables.dir = '/opt/legacy'\n"
        "config.running.mpi = 3\n",
        encoding="utf8",
    )

    migrated = migrate_legacy_configuration(legacy_path, yaml_path, definition)

    assert migrated.executables.dir() == "/opt/legacy"
    assert migrated.running.mpi() == 3
    data = load_yaml_document(yaml_path)
    assert data["profiles"]["legacy"]["executable_dir"] == "/opt/legacy"
    assert data["profiles"]["legacy"]["mpi"] == 3
    assert "mpi_warning" not in data["profiles"]["legacy"]


def test_clean_legacy_configuration_removes_top_level_assignments(tmp_path):
    path = tmp_path / "__init__.py"
    path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "\n"
        "# MPI configuration\n"
        "config.running.mpi = [\n"
        "    'mpirun',\n"
        "    '-n',\n"
        "    '4',\n"
        "]\n"
        "\n"
        "value = 42\n",
        encoding="utf8",
    )

    assert clean_legacy_configuration(path) == 1
    text = path.read_text(encoding="utf8")
    assert "config.running.mpi" not in text
    assert "# MPI configuration" not in text
    assert "from ase2sprkkr.configuration import config" not in text
    assert "value = 42" in text


def test_clean_legacy_configuration_keeps_dynamic_config_use(tmp_path):
    path = tmp_path / "__init__.py"
    path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.executables.dir = '/opt/sprkkr'\n"
        "if enabled:\n"
        "    config.running.mpi = 8\n",
        encoding="utf8",
    )

    assert clean_legacy_configuration(path) == 1
    text = path.read_text(encoding="utf8")
    assert "config.executables.dir" not in text
    assert "from ase2sprkkr.configuration import config" in text
    assert "config.running.mpi = 8" in text
