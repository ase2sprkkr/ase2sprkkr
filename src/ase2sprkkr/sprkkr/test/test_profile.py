import configparser
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr.profile import (
    SPRKKRProfile,
    profile_from_ase_command,
    resolve_sprkkr_profile,
)


def ase_config(command):
    parser = configparser.ConfigParser()
    parser["sprkkr"] = {"command": command}
    return SimpleNamespace(parser=parser)


def test_ase_command_directory_is_no_suffix():
    profile = profile_from_ase_command("/opt/sprkkr/bin")

    assert profile.executable_dir == "/opt/sprkkr/bin"
    assert profile.executable_suffix == ""
    assert profile.mpi is False


def test_ase_command_placeholder_maps_to_suffix():
    profile = profile_from_ase_command("/opt/sprkkr/bin/{}.mpi")

    assert profile.executable_dir == "/opt/sprkkr/bin"
    assert profile.executable_suffix == ".mpi"


@pytest.mark.parametrize(
    "command",
    (
        "/opt/{}/bin",
        "/opt/sprkkr/foo_{}",
        "/opt/sprkkr/{}/{}",
    ),
)
def test_invalid_ase_command_placeholder_is_rejected(command):
    with pytest.raises(ValueError, match="placeholder"):
        profile_from_ase_command(command)


def test_runtime_resolver_uses_ase_profile(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: ase_config("/opt/sprkkr/{}.x")),
    )
    configuration = definition.create_object()

    profile = resolve_sprkkr_profile(configuration=configuration)

    assert profile == SPRKKRProfile(
        executable_dir="/opt/sprkkr",
        executable_suffix=".x",
        mpi=False,
    )


def test_runtime_resolver_prefers_named_default_over_ase(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: ase_config("/ase/bin")),
    )
    configuration = definition.create_object()
    configuration.profiles.set(
        {
            "cluster": {
                "mpi": 8,
                "executable_dir": "/cluster/bin",
                "executable_suffix": ".mpi",
            }
        }
    )
    configuration.default_profile.set("cluster")

    profile = resolve_sprkkr_profile(configuration=configuration)

    assert profile == SPRKKRProfile(
        executable_dir="/cluster/bin",
        executable_suffix=".mpi",
        mpi=8,
    )


def test_runtime_resolver_accepts_explicit_profile_name(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    configuration = definition.create_object()
    configuration.profiles.set({"cluster": {"mpi": 4}})

    profile = resolve_sprkkr_profile("cluster", configuration=configuration)

    assert profile.mpi == 4
    assert profile.executable_dir is False
    assert profile.executable_suffix == ""


def test_explicit_runtime_profile_is_normalized():
    source = SPRKKRProfile(executable_dir="/bin", executable_suffix=None, mpi=None)

    profile = resolve_sprkkr_profile(source)

    assert profile == SPRKKRProfile("/bin", "", False)
    assert source == SPRKKRProfile("/bin", None, None)


def test_profile_with_overrides_returns_self_when_no_change_is_needed():
    profile = SPRKKRProfile("/cluster", ".x", 4)

    assert profile.with_overrides() is profile
    assert profile.with_overrides(mpi=None, executable_dir=None, executable_suffix=None) is profile
    assert profile.with_overrides(mpi=4, executable_dir="/cluster", executable_suffix=".x") is profile


def test_profile_with_overrides_copies_only_changed_values():
    profile = SPRKKRProfile("/cluster", ".x", 4)

    changed = profile.with_overrides(mpi=2, executable_suffix=False)

    assert changed == SPRKKRProfile("/cluster", False, 2)
    assert changed is not profile
    assert profile == SPRKKRProfile("/cluster", ".x", 4)
