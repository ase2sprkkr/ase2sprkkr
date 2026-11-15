import argparse
from functools import cache
from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr_cli.commands import config as config_command


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
        ("profiles.legacy.mpi", "8"),
        ("profiles.legacy", "{mpi: 8}"),
    ),
)
def test_config_set_rejects_reserved_legacy_profile(monkeypatch, capsys, field, value):
    config = definition.create_object()

    saved = _run_with_config(monkeypatch, config, "--set", field, value)

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
    assert "virtual profile" in capsys.readouterr().out


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

def test_ase_profile_cli_modifies_existing_config_and_not_yaml(tmp_path, monkeypatch):
    import configparser
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    path = tmp_path / "ase.ini"
    path.write_text("[espresso]\ncommand = pw.x\n\n[sprkkr]\ncommand = /old/kkrscf_v7\nmpi = 4\n")
    monkeypatch.setenv("ASE_CONFIG_PATH", str(path))
    parser = configparser.ConfigParser()
    parser.read(path)
    external = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: external),
    )
    config = definition.create_object()
    assert config.profiles.resolve_name() == "ase"

    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.default.executable_dir", "/new"
    )
    assert saved == []
    assert config.default_profile() is None
    assert config.profiles["ase"].executable_dir() == "/new"
    assert config.profiles["ase"].executable_suffix() == "_v7"
    assert "command = /new/kkrscf_v7" in path.read_text()
    assert "[espresso]\ncommand = pw.x" in path.read_text()

    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.ase.executable_suffix", "_v8"
    )
    assert saved == []
    assert config.profiles["ase"].command() == "/new/kkrscf_v8"
    assert "command = /new/kkrscf_v8" in path.read_text()


def test_ase_profile_cli_creates_config_when_missing(tmp_path, monkeypatch):
    import configparser
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    target = tmp_path / "config" / "ase.ini"
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    external = SimpleNamespace(parser=configparser.ConfigParser())
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: external),
    )
    config = definition.create_object()
    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.default.executable_dir", "/installed"
    )
    assert saved == []
    assert config.default_profile() is None
    assert "command = /installed/kkrscf" in target.read_text()
    assert config.profiles.resolve_name() == "ase"


def test_ase_profile_cli_updates_other_values_and_unsets(tmp_path, monkeypatch):
    import configparser
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    target = tmp_path / "ase.ini"
    target.write_text("[sprkkr]\ncommand = /opt/kkrscf\nmpi = 8\n")
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    parser = configparser.ConfigParser()
    parser.read(target)
    external = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: external),
    )
    config = definition.create_object()

    assert not _run_with_config(
        monkeypatch, config, "--set", "profiles.ase.mpi", "false"
    )
    assert config.profiles["ase"].mpi() is False
    assert "mpi = false" in target.read_text()
    assert not _run_with_config(
        monkeypatch, config, "--unset", "profiles.default.mpi"
    )
    assert "mpi =" not in target.read_text()
    assert config.profiles["ase"].mpi() == "auto"


def test_ase_profile_cli_invalid_value_keeps_file_and_yaml(tmp_path, monkeypatch, capsys):
    import configparser
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    target = tmp_path / "ase.ini"
    original = "[sprkkr]\ncommand = /opt/kkrscf\nmpi = 4\n"
    target.write_text(original)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    parser = configparser.ConfigParser()
    parser.read(target)
    external = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: external),
    )
    config = definition.create_object()
    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.ase.not_an_option", "4"
    )
    assert saved == []
    assert target.read_text() == original
    assert "Invalid configuration value" in capsys.readouterr().out


def test_default_profile_alias_uses_explicit_named_profile(monkeypatch):
    config = definition.create_object()
    config.profiles.set({"cluster": {"mpi": 4}}, unknown="fail", merge=True)
    config.default_profile.set("cluster")
    saved = _run_with_config(
        monkeypatch, config, "--set", "profiles.default.executable_suffix", "_v8"
    )
    assert config.default_profile() == "cluster"
    assert config.profiles["cluster"].executable_suffix() == "_v8"
    assert len(saved) == 1
    assert "default" not in saved[0]["profiles"]

def test_cli_system_flag_requires_ase_profile(monkeypatch, capsys):
    config = definition.create_object()
    saved = _run_with_config(
        monkeypatch, config, "--system",
        "--set", "profiles.cluster.mpi", "4",
    )
    assert saved == []
    assert "cluster" not in config.profiles
    assert "--system is only valid for ASE" in capsys.readouterr().out


def test_cli_system_option_is_forwarded_to_ase_writer(
    tmp_path, monkeypatch,
):
    import configparser
    from ase2sprkkr.profiles import ProfilesConfigurationContainer
    from ase2sprkkr import ase_profile_storage as storage

    system_dir = tmp_path / "etc"
    system_dir.mkdir()
    system_config = system_dir / "ase.ini"
    system_config.write_text("[sprkkr]\ncommand = /old/kkrscf\n")
    monkeypatch.setattr(storage.sys, "platform", "linux")
    monkeypatch.setattr(storage, "_SYSTEM_CONFIG_DIRECTORY", system_dir)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(system_config))
    parser = configparser.ConfigParser()
    parser.read(system_config)
    external = SimpleNamespace(parser=parser)
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: external),
    )
    config = definition.create_object()

    # Normal CLI writes cannot modify a system-only configuration.
    saved = _run_with_config(
        monkeypatch, config,
        "--set", "profiles.ase.executable_dir", "/unsafe",
    )
    assert saved == []
    assert system_config.read_text() == "[sprkkr]\ncommand = /old/kkrscf\n"

    # Explicit permission is needed for the file on this path.
    saved = _run_with_config(
        monkeypatch, config,
        "--system", "--set", "profiles.ase.executable_dir", "/new",
    )
    assert saved == []
    assert "command = /new/kkrscf" in system_config.read_text()


def test_cli_rejects_undocumented_system_wide_alias():
    with pytest.raises(SystemExit):
        _args("--system-wide", "--set", "profiles.ase.mpi", "4")


def test_cli_accepts_system_flag():
    assert _args("--system", "--set", "profiles.ase.mpi", "4").system is True


