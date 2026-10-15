import argparse
from types import SimpleNamespace

from ase2sprkkr.configuration import definition
from ase2sprkkr.tools.commands import config as config_command


def _args(*argv):
    parser = argparse.ArgumentParser()
    config_command.parser(parser)
    return parser.parse_args(argv)


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


def test_config_set_rejects_reserved_profile(monkeypatch, capsys):
    config = definition.create_object()

    saved = _run_with_config(
        monkeypatch,
        config,
        "--set",
        "profiles.ase.mpi",
        "8",
    )

    assert saved == []
    assert "reserved profile name" in capsys.readouterr().out
