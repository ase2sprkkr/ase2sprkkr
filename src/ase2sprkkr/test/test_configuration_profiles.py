import configparser
from types import SimpleNamespace

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


def ase_config(command="sprkkr"):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"command": command}
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


def test_profile_mpi_policy_has_effective_defaults_without_explicit_values():
    config = fresh_config()
    profile = config.profiles.add("local")

    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is True

    explicit = config.as_dict(only_changed="explicit", generated=False) or {}
    stored = (explicit.get("profiles") or {}).get("local", {})
    assert "mpi" not in stored
    assert "mpi_warning" not in stored


def test_named_profile_can_override_mpi_warning():
    config = fresh_config()
    profile = config.profiles.add("cluster")

    profile.mpi_warning.set(False)

    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is False
    assert config.as_dict(only_changed="explicit", generated=False)["profiles"][
        "cluster"
    ] == {"mpi_warning": False}


def test_legacy_profile_exists_only_for_explicit_legacy_configuration():
    config = fresh_config()

    assert config.running.mpi() == "auto"
    assert config.running.mpi_warning() is True
    assert "legacy" not in config.profiles

    config.running.mpi.set(8)

    assert "legacy" in config.profiles
    assert config.profiles["legacy"].mpi() == 8
    assert config.profiles["legacy"].mpi_warning() is True
    assert "legacy" in list(config.profiles)


def test_legacy_options_are_backed_by_internal_profile():
    config = fresh_config()

    config.running.mpi.set(8)
    config.running.mpi_warning.set(False)
    config.executables.dir.set("/opt/sprkkr")

    assert config.profiles.legacy_option("mpi")() == 8
    assert config.profiles.legacy_option("mpi_warning")() is False
    assert config.profiles.legacy_option("executable_dir")() == "/opt/sprkkr"
    assert config.profiles["legacy"].as_dict() == {
        "mpi": 8,
        "mpi_warning": False,
        "executable_dir": "/opt/sprkkr",
        "executable_suffix": None,
    }

    config.running.mpi.set(None)
    config.running.mpi_warning.set(None)
    config.executables.dir.set(None)
    assert "legacy" not in config.profiles
    assert "legacy" not in list(config.profiles)


def test_suffix_environment_default_does_not_create_legacy(monkeypatch):
    monkeypatch.setenv("SPRKKR_EXECUTABLE_SUFFIX", "_mpi")
    config = fresh_config()

    assert config.executables.suffix() == "_mpi"
    assert "legacy" not in config.profiles

    config.executables.suffix.set("_custom")
    assert config.profiles["legacy"].executable_suffix() == "_custom"


def test_legacy_profile_is_read_only_but_persistent():
    config = fresh_config()
    config.running.mpi.set(8)
    legacy = config.profiles["legacy"]

    with pytest.raises(TypeError, match="read-only"):
        legacy.mpi.set(4)
    with pytest.raises(TypeError, match="read-only"):
        legacy.set({"mpi": 4})
    with pytest.raises(KeyError):
        config.find("profiles.legacy.mpi")

    assert config.as_dict(only_changed="explicit")["profiles"]["legacy"]["mpi"] == 8


def test_persisted_legacy_profile_backs_deprecated_options():
    config = fresh_config()

    config.set(
        {
            "profiles": {
                "legacy": {
                    "mpi": 6,
                    "mpi_warning": False,
                    "executable_dir": "/opt/legacy",
                    "executable_suffix": "_v",
                }
            }
        },
        unknown="fail",
    )

    assert config.running.mpi() == 6
    assert config.running.mpi_warning() is False
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


def test_ase_profile_is_live_read_only_virtual_view(monkeypatch):
    external = ase_config("mpirun -n 8 /opt/sprkkr/kkrscf")
    external.parser["sprkkr"]["custom"] = "value"
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )

    config = fresh_config()

    assert "ase" in config.profiles
    profile = config.profiles["ase"]
    assert profile.command() == "mpirun -n 8 /opt/sprkkr/kkrscf"
    assert profile.custom() == "value"

    external.parser["sprkkr"]["command"] = "srun /opt/sprkkr/kkrscf"
    assert profile.command() == "srun /opt/sprkkr/kkrscf"

    with pytest.raises(TypeError, match="read-only"):
        profile.command.set("other")

    assert "ase" not in (config.as_dict().get("profiles") or {})


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
    external = ase_config("ase-command")
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )
    config = fresh_config()
    config.running.mpi.set(8)
    local = config.profiles.add("local")
    local.mpi.set(2)

    # ASE is the implicit standard and wins over the deprecated legacy config.
    assert config.profiles.resolve_name() == "ase"

    # A configured ase2sprkkr default is an explicit user choice and wins.
    config.default_profile.set("local")
    assert config.profiles.resolve() is local

    # A call-site selection is the most explicit choice of all.
    assert config.profiles.resolve("legacy").source == "legacy"


def test_unknown_explicit_default_profile_is_reported():
    config = fresh_config()
    config.default_profile.set("missing")

    with pytest.raises(KeyError, match="Unknown execution profile"):
        config.profiles.resolve()
