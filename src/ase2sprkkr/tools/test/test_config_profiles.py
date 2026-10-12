from types import SimpleNamespace

import pytest

from ase2sprkkr import configuration
from ase2sprkkr.tools.commands import config as config_command


def _args(name, value):
    return SimpleNamespace(
        path=False,
        print=False,
        show=False,
        info=False,
        set=[name, value],
        edit=False,
        default=False,
        show_default=False,
        migrate=False,
        keep_init=False,
        clean_init=False,
        delete_init=False,
        init_script=False,
        print_init_script=False,
        edit_init_script=False,
    )


def _fresh_cli(monkeypatch, tmp_path):
    config = configuration.definition.create_object()
    saves = []

    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(
        configuration, "user_preferences_file", lambda: str(tmp_path / "config.yaml")
    )
    monkeypatch.setattr(
        configuration,
        "legacy_user_preferences_file",
        lambda: str(tmp_path / "__init__.py"),
    )
    monkeypatch.setattr(
        config,
        "save_user_preferences",
        lambda: saves.append(config.as_dict(only_changed="explicit", generated=False)),
    )
    return config, saves


def test_config_set_creates_profile_from_dotted_option(monkeypatch, tmp_path):
    config, saves = _fresh_cli(monkeypatch, tmp_path)

    config_command.run(_args("profiles.cluster.mpi", "8"), None)

    assert config.profiles["cluster"].mpi() == 8
    assert len(saves) == 1


def test_config_set_creates_profile_with_mpi_warning(monkeypatch, tmp_path):
    config, saves = _fresh_cli(monkeypatch, tmp_path)

    config_command.run(_args("profiles.cluster.mpi_warning", "false"), None)

    profile = config.profiles["cluster"]
    assert profile.mpi() == "auto"
    assert profile.mpi_warning() is False
    assert saves == [{"profiles": {"cluster": {"mpi_warning": False}}}]


def test_config_set_creates_profile_from_yaml_mapping(monkeypatch, tmp_path):
    config, saves = _fresh_cli(monkeypatch, tmp_path)

    config_command.run(
        _args(
            "profiles.cluster",
            "{mpi: 8, executable_dir: /blablabla}",
        ),
        None,
    )

    profile = config.profiles["cluster"]
    assert profile.mpi() == 8
    assert profile.executable_dir() == "/blablabla"
    assert len(saves) == 1


def test_config_set_profile_mapping_updates_existing_profile(monkeypatch, tmp_path):
    config, saves = _fresh_cli(monkeypatch, tmp_path)
    profile = config.profiles.add("cluster")
    profile.executable_suffix.set(".x")

    config_command.run(
        _args(
            "profiles.cluster",
            "{mpi: 4, executable_dir: /opt/sprkkr}",
        ),
        None,
    )

    profile = config.profiles["cluster"]
    assert profile.mpi() == 4
    assert profile.executable_dir() == "/opt/sprkkr"
    assert profile.executable_suffix() == ".x"
    assert len(saves) == 1


def test_config_set_invalid_profile_is_atomic(monkeypatch, tmp_path, capsys):
    config, saves = _fresh_cli(monkeypatch, tmp_path)

    config_command.run(
        _args(
            "profiles.cluster",
            "{mpi: 8, executble_dir: /typo}",
        ),
        None,
    )

    assert "cluster" not in config.profiles
    assert saves == []
    assert "Invalid configuration value" in capsys.readouterr().out


def test_config_set_rejects_reserved_profile(monkeypatch, tmp_path, capsys):
    config, saves = _fresh_cli(monkeypatch, tmp_path)

    config_command.run(
        _args("profiles.ase", "{mpi: 8}"),
        None,
    )

    assert "ase" not in config.profiles._values
    assert saves == []
    assert "reserved profile name" in capsys.readouterr().out
