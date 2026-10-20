import configparser
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer


def fresh_config():
    return definition.create_object()


def _with_ase_profile(monkeypatch):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"command": "kkrscf"}
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )


def test_profile_collection_exposes_only_public_views(monkeypatch):
    _with_ase_profile(monkeypatch)
    config = fresh_config()
    local = config.profiles.add("local")
    local.mpi.set(2)
    config.running.mpi.set(4)

    names = list(config.profiles)
    items = dict(config.profiles.items())
    values = list(config.profiles.values())

    assert names == ["local", "legacy", "ase"]
    assert list(items) == names
    assert len(config.profiles) == len(names)
    assert values[0] is local
    assert local.read_only is False
    assert config.profiles["local"].read_only is False
    assert [profile.source for profile in values[1:]] == ["legacy", "ase"]

    assert items["local"] is local
    assert items["legacy"].read_only is True
    assert items["ase"].read_only is True

    with pytest.raises(TypeError, match="read-only"):
        items["legacy"].mpi.set(8)


def test_reserved_profiles_cannot_be_created_publicly():
    config = fresh_config()

    with pytest.raises(ValueError, match="reserved profile"):
        config.profiles.add("legacy")
    with pytest.raises(ValueError, match="reserved profile"):
        config.profiles.add("ase")


def test_reserved_profiles_cannot_be_set_publicly():
    config = fresh_config()

    for name in ("legacy", "ase"):
        with pytest.raises(ValueError, match="read-only"):
            config.profiles.set({name: {"mpi": 8}}, unknown="fail", merge=True)
        with pytest.raises(ValueError, match="read-only"):
            config.set(
                {"profiles": {name: {"mpi": 8}}},
                unknown="fail",
            )


def test_public_profile_replace_preserves_legacy_backing():
    config = fresh_config()
    config.running.mpi.set(4)
    config.profiles.add("old").mpi.set(2)

    config.profiles.set({"new": {"mpi": 3}}, unknown="fail")

    assert list(config.profiles) == ["new", "legacy"]
    assert config.profiles["new"].mpi() == 3
    assert config.profiles["legacy"].mpi() == 4


def test_generic_lookup_does_not_expose_legacy_backing():
    config = fresh_config()
    config.running.mpi.set(4)

    with pytest.raises(KeyError):
        config.profiles.get_member("legacy.mpi", unknown="fail")
