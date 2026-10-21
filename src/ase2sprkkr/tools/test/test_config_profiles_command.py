import argparse
from functools import cache
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.tools.commands import config as config_command


@cache
def _parser():
    parser = argparse.ArgumentParser()
    config_command.parser(parser)
    return parser


def _args(*argv):
    return _parser().parse_args(argv)


def _run_with_config(monkeypatch, config, *argv):
    import ase2sprkkr.configuration as configuration

    saved = []
    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(
        type(config),
        "save_user_preferences",
        lambda self: saved.append(
            self.as_dict(only_changed="explicit", generated=False)
        ),
    )
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: "/tmp/unused-config.yaml")
    monkeypatch.setattr(configuration, "legacy_user_preferences_file", lambda: "/tmp/unused-init.py")

    config_command.run(_args(*argv), SimpleNamespace())
    return saved


def test_config_set_creates_profile_field(monkeypatch):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch,
        config,
        "--set",
        "profiles.cluster.mpi",
        "8",
    )

    assert config.profiles["cluster"].mpi() == 8
    assert len(saved) == 1


def test_config_set_creates_whole_profile_from_yaml_mapping(monkeypatch):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch,
        config,
        "--set",
        "profiles.cluster",
        "{mpi: 8, executable_dir: /blablabla}",
    )

    profile = config.profiles["cluster"]
    assert profile.mpi() == 8
    assert profile.executable_dir() == "/blablabla"
    assert len(saved) == 1


def test_config_set_profile_is_atomic_on_invalid_member(monkeypatch, capsys):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch,
        config,
        "--set",
        "profiles.cluster",
        "{mpi: 8, executble_dir: /typo}",
    )

    assert "cluster" not in config.profiles
    assert saved == []
    assert "Invalid configuration value" in capsys.readouterr().out


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("profiles.ase.mpi", "8"),
        ("profiles.ase", "{mpi: 8}"),
    ),
)
def test_config_set_rejects_reserved_profile(monkeypatch, capsys, field, value):
    config = definition.create_object()

    saved = _run_with_config(monkeypatch, config, "--set", field, value)

    assert "ase" not in config.profiles._values
    assert saved == []
    assert "reserved profile name" in capsys.readouterr().out


def test_config_set_creates_profile_with_mpi_warning(monkeypatch):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.cluster.mpi_warning", "false",
    )

    profile = config.profiles["cluster"]
    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is False
    assert saved == [{"profiles": {"cluster": {"mpi_warning": False}}}]


def test_config_set_profile_mapping_updates_existing_profile(monkeypatch):
    config = definition.create_object()
    profile = config.profiles.add("cluster")
    profile.executable_suffix.set(".x")

    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.cluster",
        "{mpi: 4, executable_dir: /opt/sprkkr}",
    )

    profile = config.profiles["cluster"]
    assert profile.mpi() == 4
    assert profile.executable_dir() == "/opt/sprkkr"
    assert profile.executable_suffix() == ".x"
    assert len(saved) == 1


def test_config_unset_removes_named_profile(monkeypatch):
    config = definition.create_object()
    config.profiles.set(
        {
            "cluster": {"mpi": 8},
            "local": {"mpi": 4},
        },
        unknown="fail",
    )

    saved = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "profiles.cluster",
    )

    assert "cluster" not in config.profiles
    assert config.profiles["local"].mpi() == 4
    assert saved == [{"profiles": {"local": {"mpi": 4}}}]


def test_config_unset_profile_field_preserves_empty_default_profile(monkeypatch):
    config = definition.create_object()
    config.profiles.set({"cluster": {"mpi": 8}}, unknown="fail")
    config.default_profile.set("cluster")

    saved = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "profiles.cluster.mpi",
    )

    assert config.profiles["cluster"].mpi() == "auto"
    assert saved == [{"default_profile": "cluster", "profiles": {"cluster": {}}}]


def test_config_unset_rejects_reserved_profile(monkeypatch, capsys):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "profiles.ase",
    )

    assert saved == []
    assert "reserved profile name" in capsys.readouterr().out


def test_config_unset_rejects_explicit_default_profile(monkeypatch, capsys):
    config = definition.create_object()
    config.profiles.set({"cluster": {"mpi": 8}}, unknown="fail")
    config.default_profile.set("cluster")

    saved = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "profiles.cluster",
    )

    assert "cluster" in config.profiles
    assert config.default_profile() == "cluster"
    assert saved == []
    assert "default execution profile" in capsys.readouterr().out


def test_config_unset_default_profile_allows_later_profile_removal(monkeypatch):
    config = definition.create_object()
    config.profiles.set({"cluster": {"mpi": 8}}, unknown="fail")
    config.default_profile.set("cluster")

    first_save = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "default_profile",
    )
    second_save = _run_with_config(
        monkeypatch,
        config,
        "--unset",
        "profiles.cluster",
    )

    assert config.default_profile() is None
    assert "cluster" not in config.profiles
    assert len(first_save) == 1
    assert len(second_save) == 1


def test_config_list_profiles_marks_explicit_default(monkeypatch, capsys):
    config = definition.create_object()
    config.profiles.set(
        {
            "cluster": {"mpi": 8},
            "local": {"mpi": 4},
        },
        unknown="fail",
    )
    config.default_profile.set("cluster")

    saved = _run_with_config(monkeypatch, config, "--list-profiles")

    output = capsys.readouterr().out.splitlines()
    assert "cluster (default)" in output
    assert "local" in output
    assert saved == []
