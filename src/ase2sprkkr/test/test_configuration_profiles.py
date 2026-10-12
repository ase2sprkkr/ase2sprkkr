import configparser
from types import SimpleNamespace
import warnings

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer


@pytest.fixture(autouse=True)
def no_host_ase_profile(monkeypatch):
    """Keep profile tests independent of the developer machine's ASE config."""
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )


def fresh_config():
    return definition.create_object()


def ase_config(command="kkrscf", **values):
    parser = configparser.ConfigParser()
    section = {"command": command}
    section.update({name: str(value) for name, value in values.items()})
    parser["sprkkr"] = section
    return SimpleNamespace(parser=parser)


def test_named_profile_is_persistent_configuration():
    config = fresh_config()

    profile = config.profiles.add("local")
    profile.mpi.set(4)
    profile.executable_dir.set("/opt/sprkkr/bin")

    assert list(config.profiles) == ["local"]
    assert config.profiles["local"] is profile
    assert profile.mpi() == 4
    assert config.as_dict()["profiles"]["local"]["mpi"] == 4


def test_legacy_profile_exists_only_for_explicit_legacy_configuration():
    config = fresh_config()

    assert "legacy" not in config.profiles

    with pytest.warns(FutureWarning):
        config.running.mpi.set(8)

    assert "legacy" in config.profiles
    assert config.profiles["legacy"].mpi() == 8
    assert "legacy" in list(config.profiles)


def test_legacy_options_are_backed_by_internal_profile():
    config = fresh_config()

    with pytest.warns(FutureWarning):
        config.running.mpi.set(8)
    with pytest.warns(FutureWarning):
        config.executables.dir.set("/opt/sprkkr")

    assert config.profiles.legacy_option("mpi")() == 8
    assert config.profiles.legacy_option("executable_dir")() == "/opt/sprkkr"
    assert config.profiles["legacy"].as_dict() == {
        "mpi": 8,
        "executable_dir": "/opt/sprkkr",
        "executable_suffix": None,
    }

    with pytest.warns(FutureWarning):
        config.running.mpi.set(None)
    with pytest.warns(FutureWarning):
        config.executables.dir.set(None)
    assert "legacy" not in config.profiles
    assert "legacy" not in list(config.profiles)


def test_suffix_environment_default_does_not_create_legacy(monkeypatch):
    monkeypatch.setenv("SPRKKR_EXECUTABLE_SUFFIX", "_mpi")
    config = fresh_config()

    assert config.executables.suffix() == "_mpi"
    assert "legacy" not in config.profiles

    with pytest.warns(FutureWarning):
        config.executables.suffix.set("_custom")
    assert config.profiles["legacy"].executable_suffix() == "_custom"


def test_legacy_profile_is_read_only_but_persistent():
    config = fresh_config()
    with pytest.warns(FutureWarning):
        config.running.mpi.set(8)
    legacy = config.profiles["legacy"]

    with pytest.raises(TypeError, match="read-only"):
        legacy.mpi.set(4)
    with pytest.raises(TypeError, match="read-only"):
        legacy.set({"mpi": 4})
    with pytest.raises(KeyError):
        config.find("profiles.legacy.mpi")

    assert config.as_dict(only_changed="explicit")["profiles"]["legacy"]["mpi"] == 8


def test_legacy_write_warning_shows_both_profile_migrations():
    config = fresh_config()

    with pytest.warns(FutureWarning) as caught:
        config.executables.dir.set("/opt/sprkkr")

    message = str(caught[0].message)
    assert "config.executables.dir is deprecated" in message
    assert "[sprkkr]" in message
    assert "command = /opt/sprkkr/kkrscf_v8" in message
    assert "mpi = 8" in message
    assert "ase2sprkkr config --set profiles.local" in message
    assert "ase2sprkkr config --set default_profile local" in message


def test_legacy_write_warning_is_emitted_once():
    config = fresh_config()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default", FutureWarning)
        config.running.mpi.set(8)
        config.executables.dir.set("/opt/sprkkr")
        config.executables.suffix.set("_v8")

    assert len(caught) == 1


def test_persisted_legacy_profile_backs_deprecated_options():
    config = fresh_config()

    config.set(
        {
            "profiles": {
                "legacy": {
                    "mpi": 6,
                    "executable_dir": "/opt/legacy",
                    "executable_suffix": "_v",
                }
            }
        },
        unknown="fail",
    )

    assert config.running.mpi() == 6
    assert config.executables.dir() == "/opt/legacy"
    assert config.executables.suffix() == "_v"
    with pytest.raises(TypeError, match="read-only"):
        config.profiles["legacy"].mpi.set(2)


def test_reserved_virtual_profile_names_cannot_be_added():
    config = fresh_config()

    with pytest.raises(ValueError, match="reserved profile name"):
        config.profiles.add("legacy")
    with pytest.raises(ValueError, match="reserved profile name"):
        config.profiles.add("ase")


def test_ase_profile_uses_common_profile_schema(monkeypatch):
    external = ase_config(
        "/opt/sprkkr/kkrscf.mpi",
        mpi=8,
        executor="ssh",
        host="user@cluster",
        input_sync="always",
    )
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )

    config = fresh_config()
    profile = config.profiles["ase"]

    assert profile.command() == "/opt/sprkkr/kkrscf.mpi"
    assert profile.executable_dir() == "/opt/sprkkr"
    assert profile.executable_suffix() == ".mpi"
    assert profile.mpi() == 8
    assert profile.executor() == "ssh"
    assert profile.host() == "user@cluster"
    assert profile.input_sync() == "always"

    external.parser["sprkkr"]["mpi"] = "4"
    assert profile.mpi() == 4

    with pytest.raises(TypeError, match="read-only"):
        profile.mpi.set(2)

    assert "ase" not in (config.as_dict().get("profiles") or {})


@pytest.mark.parametrize(
    "command",
    (
        "mpirun -n 8 /opt/sprkkr/kkrscf",
        "srun /opt/sprkkr/kkrscf",
        "/opt/sprkkr/kkrspec",
        "/opt/sprkkr/{}.mpi",
    ),
)
def test_ase_profile_rejects_noncanonical_command(monkeypatch, command):
    external = ase_config(command)
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )
    config = fresh_config()

    with pytest.raises(ValueError, match="kkrscf|prefixes|wrappers"):
        config.profiles["ase"].configuration_profile()


def test_ase_profile_reuses_common_validation(monkeypatch):
    external = ase_config("/opt/sprkkr/kkrscf", executor="not-an-executor")
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )
    config = fresh_config()

    with pytest.raises(ValueError):
        config.profiles["ase"].executor()


def test_incomplete_ase_section_is_not_an_available_profile(monkeypatch):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"other": "value"}
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )

    config = fresh_config()

    assert "ase" not in config.profiles
    with pytest.raises(KeyError):
        config.profiles["ase"]


def test_profile_resolution_precedence(monkeypatch):
    external = ase_config("kkrscf")
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )
    config = fresh_config()
    with pytest.warns(FutureWarning):
        config.running.mpi.set(8)
    local = config.profiles.add("local")
    local.mpi.set(2)

    assert config.profiles.resolve_name() == "ase"

    config.default_profile.set("local")
    assert config.profiles.resolve() is local

    assert config.profiles.resolve("legacy").source == "legacy"


def test_unknown_explicit_default_profile_is_reported():
    config = fresh_config()
    config.default_profile.set("missing")

    with pytest.raises(KeyError, match="Unknown execution profile"):
        config.profiles.resolve()
