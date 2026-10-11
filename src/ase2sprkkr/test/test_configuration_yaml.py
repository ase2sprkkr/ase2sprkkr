from pathlib import Path

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_yaml import (
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
    profile.executable_dir.set("/opt/sprkkr")
    config.default_profile.set("local")

    save_configuration_yaml(config, path)

    loaded = fresh_config()
    load_configuration_yaml(loaded, path)

    assert loaded.default_profile() == "local"
    assert loaded.profiles["local"].mpi() == 4
    assert loaded.profiles["local"].executable_dir() == "/opt/sprkkr"


def test_yaml_does_not_persist_virtual_legacy_profile(tmp_path):
    path = tmp_path / "config.yaml"
    config = fresh_config()
    config.running.mpi.set(8)

    save_configuration_yaml(config, path)
    data = load_yaml_document(path)

    assert "legacy" not in (data.get("profiles") or {})
    assert "running" not in data or "mpi" not in data["running"]


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


def test_yaml_rejects_unknown_configuration_keys(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("unknown_option: value\n", encoding="utf8")

    with pytest.raises((KeyError, ValueError)):
        load_configuration_yaml(fresh_config(), path)
