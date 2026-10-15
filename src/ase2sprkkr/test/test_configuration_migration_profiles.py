import configparser
from types import SimpleNamespace
import warnings

from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_migration import migrate_legacy_configuration
from ase2sprkkr.configuration_yaml import load_yaml_document
from ase2sprkkr.profiles import ProfilesConfigurationContainer


def _write_legacy(path, body):
    path.write_text(
        "from ase2sprkkr.configuration import config\n" + body,
        encoding="utf8",
    )


def test_migration_converts_legacy_execution_settings_to_default_profile(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    source = tmp_path / "__init__.py"
    destination = tmp_path / "config.yaml"
    _write_legacy(
        source,
        "config.executables.dir = '/opt/sprkkr'\n"
        "config.executables.suffix = '_v8'\n"
        "config.running.mpi = 4\n",
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", FutureWarning)
        migrated = migrate_legacy_configuration(source, destination, definition)

    assert caught == []
    assert migrated.default_profile() == "local"
    assert "legacy" not in migrated.profiles
    profile = migrated.profiles["local"]
    assert profile.mpi() == 4
    assert profile.executable_dir() == "/opt/sprkkr"
    assert profile.executable_suffix() == "_v8"

    data = load_yaml_document(destination)
    assert data["default_profile"] == "local"
    assert "legacy" not in data["profiles"]
    assert data["profiles"]["local"] == {
        "mpi": 4,
        "executable_dir": "/opt/sprkkr",
        "executable_suffix": "_v8",
    }


def test_migration_keeps_ase_profile_precedence(tmp_path, monkeypatch):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"command": "/opt/ase/kkrscf", "mpi": "2"}
    ase_config = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: ase_config),
    )
    source = tmp_path / "__init__.py"
    destination = tmp_path / "config.yaml"
    _write_legacy(
        source,
        "config.executables.dir = '/opt/legacy'\n"
        "config.running.mpi = 8\n",
    )

    migrated = migrate_legacy_configuration(source, destination, definition)

    assert migrated.default_profile() is None
    assert migrated.profiles.resolve_name() == "ase"
    assert "legacy" not in migrated.profiles
    assert migrated.profiles["local"].mpi() == 8
    assert migrated.profiles["local"].executable_dir() == "/opt/legacy"

    data = load_yaml_document(destination)
    assert "default_profile" not in data
    assert "legacy" not in data["profiles"]
    assert data["profiles"]["local"]["mpi"] == 8


def test_migration_does_not_overwrite_existing_local_profile(tmp_path, monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    source = tmp_path / "__init__.py"
    destination = tmp_path / "config.yaml"
    _write_legacy(
        source,
        "config.profiles.add('local').mpi.set(2)\n"
        "config.executables.dir = '/opt/legacy'\n"
        "config.running.mpi = 8\n",
    )

    migrated = migrate_legacy_configuration(source, destination, definition)

    assert migrated.profiles["local"].mpi() == 2
    assert migrated.profiles["local-2"].mpi() == 8
    assert migrated.profiles["local-2"].executable_dir() == "/opt/legacy"
    assert migrated.default_profile() == "local-2"
