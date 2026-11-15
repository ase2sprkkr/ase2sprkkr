"""ASE config.ini selection and atomic updates, independent of SPR-KKR binaries."""

import os
from pathlib import Path
import stat

import pytest

from ase2sprkkr import ase_profile_storage as storage


def _system_and_user_paths(tmp_path, monkeypatch):
    """Simulate Linux system and user config roots without writing to /etc."""
    system_dir = tmp_path / "etc"
    user_dir = tmp_path / "home" / "user"
    system_dir.mkdir(parents=True)
    user_dir.mkdir(parents=True)
    monkeypatch.setattr(storage.sys, "platform", "linux")
    monkeypatch.setattr(storage, "_SYSTEM_CONFIG_DIRECTORY", system_dir)
    system_file = system_dir / "ase.ini"
    user_file = user_dir / "ase.ini"
    monkeypatch.setenv("ASE_CONFIG_PATH", f"{system_file}:{user_file}")
    return system_file, user_file


def test_linux_system_sprkkr_is_not_modified_by_default(tmp_path, monkeypatch):
    system_file, user_file = _system_and_user_paths(tmp_path, monkeypatch)
    original = "[sprkkr]\ncommand = /system/kkrscf\n"
    system_file.write_text(original)
    assert storage.ase_config_target() == user_file

    storage.write_ase_sprkkr({"command": "/user/kkrscf"})
    assert system_file.read_text() == original
    assert user_file.read_text() == "[sprkkr]\ncommand = /user/kkrscf\n"


def test_linux_system_config_opt_in(tmp_path, monkeypatch):
    system_file, user_file = _system_and_user_paths(tmp_path, monkeypatch)
    system_file.write_text("[sprkkr]\ncommand = /system/kkrscf\n")

    assert storage.ase_config_target(system=True) == system_file
    storage.write_ase_sprkkr({"command": "/updated/kkrscf"}, system=True)
    assert "command = /updated/kkrscf" in system_file.read_text()
    assert not user_file.exists()


def test_linux_only_system_config_path_raises(tmp_path, monkeypatch):
    system_file, user_file = _system_and_user_paths(tmp_path, monkeypatch)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(system_file))
    original = "[sprkkr]\ncommand = /system/kkrscf\n"
    system_file.write_text(original)
    with pytest.raises(PermissionError, match="system=True"):
        storage.ase_config_target()
    with pytest.raises(PermissionError, match="system=True"):
        storage.write_ase_sprkkr({"command": "/unsafe/kkrscf"})
    assert system_file.read_text() == original
    assert not user_file.exists()


def test_linux_user_symlink_into_system_config_is_rejected(tmp_path, monkeypatch):
    system_file, user_file = _system_and_user_paths(tmp_path, monkeypatch)
    original = "[sprkkr]\ncommand = /system/kkrscf\n"
    system_file.write_text(original)
    user_file.symlink_to(system_file)
    with pytest.raises(PermissionError, match="system=True"):
        storage.ase_config_target()
    assert system_file.read_text() == original


def test_linux_higher_priority_system_config_requires_opt_in(tmp_path, monkeypatch):
    system_file, user_file = _system_and_user_paths(tmp_path, monkeypatch)
    monkeypatch.setenv("ASE_CONFIG_PATH", f"{user_file}:{system_file}")
    user_file.write_text("[sprkkr]\ncommand = /user/kkrscf\n")
    system_file.write_text("[sprkkr]\ncommand = /system/kkrscf\n")

    # Editing the user file would not take effect, because the system
    # section wins when ASE merges the configured paths.
    with pytest.raises(PermissionError, match="higher-priority system"):
        storage.ase_config_target()
    assert storage.ase_config_target(system=True) == system_file


def test_windows_programdata_is_system_config_but_user_home_is_not(monkeypatch):
    monkeypatch.setattr(storage.sys, "platform", "win32")
    assert storage._is_system_config_path(Path(r"C:\\ProgramData\\ASE\\config.ini"))
    assert storage._is_system_config_path(Path(r"C:\\Users\\Default\\config.ini"))
    assert not storage._is_system_config_path(Path(r"C:\\Users\\Alice\\.config\\ase.ini"))


def test_staged_ase_option_is_not_saved_until_explicit_persistence(tmp_path, monkeypatch):
    import configparser
    from types import SimpleNamespace
    from ase2sprkkr.configuration import definition
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    target = tmp_path / "ase.ini"
    original = "[sprkkr]\ncommand = /old/kkrscf\n"
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
    config.profiles.ase.executable_dir.set("/new")
    config.profiles.ase.mpi.set(4)
    assert config.profiles.ase.executable_dir() == "/new"
    assert config.profiles.ase.mpi() == 4
    assert target.read_text() == original

    config.profiles.ase.mpi.store_permanent()
    assert "command = /new/kkrscf" in target.read_text()
    assert "mpi = 4" in target.read_text()

    config.profiles.ase.mpi.set_permanent(8)
    assert "mpi = 8" in target.read_text()


def test_root_save_user_preferences_persists_staged_ase_values(tmp_path, monkeypatch):
    import configparser
    from types import SimpleNamespace
    import ase2sprkkr.configuration as configuration
    import ase2sprkkr.configuration_yaml as yaml
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    target = tmp_path / "ase.ini"
    original = "[sprkkr]\ncommand = /old/kkrscf\n"
    target.write_text(original)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    parser = configparser.ConfigParser()
    parser.read(target)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )
    config = configuration.definition.create_object()
    calls = []
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(tmp_path / "empty.yaml"))
    monkeypatch.setattr(yaml, "save_configuration_yaml", lambda *_: calls.append(True))
    config.profiles.ase.executable_suffix.set("_v9")
    assert target.read_text() == original

    config.save_user_preferences()
    assert "command = /old/kkrscf_v9" in target.read_text()
    assert calls == [True]


def test_root_save_system_flag_is_forwarded_after_in_memory_edit(
    tmp_path, monkeypatch
):
    import configparser
    from types import SimpleNamespace
    import ase2sprkkr.configuration as configuration
    import ase2sprkkr.configuration_yaml as yaml
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    system_file, _user_file = _system_and_user_paths(tmp_path, monkeypatch)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(system_file))
    original = "[sprkkr]\ncommand = /old/kkrscf\n"
    system_file.write_text(original)

    parser = configparser.ConfigParser()
    parser.read(system_file)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )
    config = configuration.definition.create_object()
    monkeypatch.setattr(
        configuration, "user_preferences_file",
        lambda: str(tmp_path / "user.yaml"),
    )
    calls = []
    monkeypatch.setattr(
        yaml, "save_configuration_yaml", lambda *_: calls.append(True),
    )

    config.profiles.ase.executable_dir.set("/new")
    assert system_file.read_text() == original
    with pytest.raises(PermissionError, match="system=True"):
        config.save_user_preferences()
    assert system_file.read_text() == original
    assert calls == []

    config.save_user_preferences(system=True)
    assert "command = /new/kkrscf" in system_file.read_text()
    assert calls == [True]


def test_set_permanent_system_flag_is_forwarded(tmp_path, monkeypatch):
    import configparser
    from types import SimpleNamespace
    from ase2sprkkr.configuration import definition
    from ase2sprkkr.profiles import ProfilesConfigurationContainer

    system_file, _user_file = _system_and_user_paths(tmp_path, monkeypatch)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(system_file))
    system_file.write_text("[sprkkr]\ncommand = /old/kkrscf\n")
    parser = configparser.ConfigParser()
    parser.read(system_file)
    monkeypatch.setattr(
        ProfilesConfigurationContainer, "_ase_config",
        staticmethod(lambda: SimpleNamespace(parser=parser)),
    )
    config = definition.create_object()
    config.profiles.ase.mpi.set_permanent(4, system=True)
    assert "mpi = 4" in system_file.read_text()


def test_edit_existing_sprkkr_profile_in_precedence_order(tmp_path, monkeypatch):
    first, second, third = (
        tmp_path / name for name in ("first.ini", "second.ini", "third.ini")
    )
    first.write_text("[sprkkr]\ncommand = /first/kkrscf\n")
    second.write_text("[espresso]\ncommand = pw.x\n")
    third.write_text("[other]\nsetting = value\n")
    monkeypatch.setenv(
        "ASE_CONFIG_PATH", os.pathsep.join(map(str, (first, second, third)))
    )

    # The highest-priority file exists but does not define [sprkkr].
    assert storage.ase_config_target() == first
    storage.write_ase_sprkkr({"command": "/updated/kkrscf"})
    assert "command = /updated/kkrscf" in first.read_text()
    assert second.read_text() == "[espresso]\ncommand = pw.x\n"
    assert third.read_text() == "[other]\nsetting = value\n"

    # Another [sprkkr] section at a higher priority wins, even when it
    # supplies only one setting and inherits command from the first file.
    second.write_text("[sprkkr]\nmpi = 4\n")
    assert storage.ase_config_target() == second
    storage.write_ase_sprkkr({"command": "/second/kkrscf"})
    assert "command = /second/kkrscf" in second.read_text()
    assert "command = /updated/kkrscf" in first.read_text()
    assert third.read_text() == "[other]\nsetting = value\n"


def test_no_sprkkr_section_uses_highest_priority_path_even_if_others_exist(
    tmp_path, monkeypatch
):
    first, second, third = (
        tmp_path / name for name in ("first.ini", "second.ini", "third.ini")
    )
    first.write_text("[espresso]\ncommand = pw.x\n")
    second.write_text("[abinit]\ncommand = abinit\n")
    monkeypatch.setenv(
        "ASE_CONFIG_PATH", os.pathsep.join(map(str, (first, second, third)))
    )

    assert storage.ase_config_target() == third
    storage.write_ase_sprkkr({"command": "/new/kkrscf"})
    assert third.read_text() == "[sprkkr]\ncommand = /new/kkrscf\n"
    assert first.read_text() == "[espresso]\ncommand = pw.x\n"
    assert second.read_text() == "[abinit]\ncommand = abinit\n"


def test_create_highest_priority_only_when_none_exists(tmp_path, monkeypatch):
    first, second, last = (tmp_path / name for name in ("a.ini", "b.ini", "c.ini"))
    monkeypatch.setenv("ASE_CONFIG_PATH", os.pathsep.join(map(str, (first, second, last))))
    assert storage.ase_config_target() == last
    storage.write_ase_sprkkr({"command": "/sprkkr/kkrscf"})
    assert last.is_file()
    assert not first.exists() and not second.exists()
    assert "command = /sprkkr/kkrscf" in last.read_text()


def test_no_configured_path_does_not_write_somewhere_else(monkeypatch):
    monkeypatch.setenv("ASE_CONFIG_PATH", "")
    with pytest.raises(ValueError, match="empty"):
        storage.ase_config_target()


def test_existing_ini_sections_comments_and_permissions_preserved(tmp_path, monkeypatch):
    original = (
        "# my ASE settings\n[espresso]\ncommand = pw.x\n"
        "\n[sprkkr]\n# user comment\ncommand = /old/kkrscf\nmpi = auto\n"
        "\n[abinit]\ncommand = abinit\n"
    )
    target = tmp_path / "ase.ini"
    target.write_text(original)
    target.chmod(0o640)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    storage.write_ase_sprkkr({"command": "/new/kkrscf_v8", "mpi": "4"})
    data = target.read_text()
    assert "[espresso]\ncommand = pw.x" in data
    assert "[abinit]\ncommand = abinit" in data
    assert "# user comment" in data
    assert "command = /new/kkrscf_v8" in data
    assert "mpi = 4" in data
    assert stat.S_IMODE(target.stat().st_mode) == 0o640


def test_symlink_is_preserved_on_atomic_write(tmp_path, monkeypatch):
    target = tmp_path / "real.ini"
    target.write_text("[sprkkr]\ncommand = kkrscf\n")
    link = tmp_path / "config.ini"
    link.symlink_to(target)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(link))
    storage.write_ase_sprkkr({"command": "/path/kkrscf"})
    assert link.is_symlink()
    assert target.read_text().strip().endswith("command = /path/kkrscf")


def test_failed_replace_preserves_original(tmp_path, monkeypatch):
    target = tmp_path / "ase.ini"
    original = "[sprkkr]\ncommand = kkrscf\n"
    target.write_text(original)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))

    def fail_replace(*args):
        raise OSError("simulated rename failure")

    monkeypatch.setattr(storage.os, "replace", fail_replace)
    with pytest.raises(OSError, match="simulated"):
        storage.write_ase_sprkkr({"command": "kkrscf_v9"})
    assert target.read_text() == original
    assert list(tmp_path.iterdir()) == [target]


def test_invalid_existing_ini_is_not_overwritten(tmp_path, monkeypatch):
    target = tmp_path / "ase.ini"
    original = "[sprkkr\ncommand = invalid\n"
    target.write_text(original)
    monkeypatch.setenv("ASE_CONFIG_PATH", str(target))
    with pytest.raises(Exception):
        storage.write_ase_sprkkr({"command": "kkrscf"})
    assert target.read_text() == original


def test_append_missing_option_after_unterminated_line():
    original = "[sprkkr]\ncommand = kkrscf"
    updated = storage._rewrite_section(original, {"mpi": "4"})
    assert "command = kkrscf\nmpi = 4" in updated

def test_multiple_ase_files_merge_and_later_values_win(tmp_path, monkeypatch):
    import configparser
    earlier = tmp_path / "first.ini"
    later = tmp_path / "second.ini"
    earlier.write_text(
        "[sprkkr]\ncommand = /first/kkrscf\nmpi = 2\n"
        "[espresso]\ncommand = pw.x\n"
    )
    later.write_text("[sprkkr]\nmpi = 4\n")
    monkeypatch.setenv("ASE_CONFIG_PATH", str(earlier) + ":" + str(later))

    parser = configparser.ConfigParser()
    parser.read([earlier, later])
    assert parser["sprkkr"]["command"] == "/first/kkrscf"
    assert parser["sprkkr"]["mpi"] == "4"
    assert storage.ase_config_target() == later
    storage.write_ase_sprkkr({"command": "/updated/kkrscf"})
    parser = configparser.ConfigParser()
    parser.read([earlier, later])
    assert parser["sprkkr"]["command"] == "/updated/kkrscf"
    assert parser["sprkkr"]["mpi"] == "4"
    assert "[espresso]" in earlier.read_text()
