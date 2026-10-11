import configparser
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer


def fresh_config():
    return definition.create_object()


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

    config.running.mpi.set(8)

    assert "legacy" in config.profiles
    assert config.profiles["legacy"].mpi() == 8
    assert "legacy" in list(config.profiles)


def test_legacy_profile_is_read_only_and_not_serialized():
    config = fresh_config()
    config.running.mpi.set(8)
    legacy = config.profiles["legacy"]

    with pytest.raises(TypeError, match="read-only"):
        legacy.mpi.set(4)
    with pytest.raises(TypeError, match="read-only"):
        legacy.set({"mpi": 4})

    profiles = config.as_dict().get("profiles") or {}
    assert "legacy" not in profiles


def test_reserved_virtual_profile_names_cannot_be_added():
    config = fresh_config()

    with pytest.raises(ValueError, match="reserved profile name"):
        config.profiles.add("legacy")
    with pytest.raises(ValueError, match="reserved profile name"):
        config.profiles.add("ase")


def test_ase_profile_is_live_read_only_virtual_view(monkeypatch):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {
        "command": "mpirun -n 8 /opt/sprkkr/kkrscf",
        "custom": "value",
    }
    ase_config = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: ase_config),
    )

    config = fresh_config()

    assert "ase" in config.profiles
    profile = config.profiles["ase"]
    assert profile.command() == "mpirun -n 8 /opt/sprkkr/kkrscf"
    assert profile.custom() == "value"

    parser["sprkkr"]["command"] = "srun /opt/sprkkr/kkrscf"
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
