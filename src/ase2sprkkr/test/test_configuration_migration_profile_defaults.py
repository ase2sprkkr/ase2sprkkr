from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_migration import migrate_legacy_configuration
from ase2sprkkr.configuration_yaml import load_yaml_document
from ase2sprkkr.profiles import ProfilesConfigurationContainer


def _write_legacy(path, body):
    path.write_text(
        "from ase2sprkkr.configuration import config\n" + body,
        encoding="utf8",
    )


def _disable_ase_profile(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )


def test_migration_does_not_freeze_mpi_defaults(tmp_path, monkeypatch):
    _disable_ase_profile(monkeypatch)
    source = tmp_path / "__init__.py"
    destination = tmp_path / "config.yaml"
    _write_legacy(
        source,
        "config.executables.dir = '/opt/sprkkr'\n"
        "config.running.mpi = 'auto'\n"
        "config.running.mpi_warning = True\n",
    )

    migrated = migrate_legacy_configuration(source, destination, definition)

    profile = migrated.profiles["local"]
    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is True
    assert load_yaml_document(destination)["profiles"]["local"] == {
        "executable_dir": "/opt/sprkkr",
    }


def test_migration_keeps_nondefault_mpi_warning(tmp_path, monkeypatch):
    _disable_ase_profile(monkeypatch)
    source = tmp_path / "__init__.py"
    destination = tmp_path / "config.yaml"
    _write_legacy(
        source,
        "config.executables.dir = '/opt/sprkkr'\n"
        "config.running.mpi_warning = False\n",
    )

    migrate_legacy_configuration(source, destination, definition)

    assert load_yaml_document(destination)["profiles"]["local"] == {
        "mpi_warning": False,
        "executable_dir": "/opt/sprkkr",
    }
