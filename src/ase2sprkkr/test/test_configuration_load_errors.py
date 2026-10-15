import pytest

import ase2sprkkr.configuration as configuration
import ase2sprkkr.configuration_yaml as configuration_yaml
from ase2sprkkr.configuration import definition


def fresh_config():
    return definition.create_object()


def test_user_preferences_propagate_yaml_loader_programming_error(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"
    yaml_path.write_text("profiles: {}\n", encoding="utf8")

    monkeypatch.setattr(configuration, "config", fresh_config())
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(
        configuration,
        "legacy_user_preferences_file",
        lambda: str(legacy_path),
    )

    def broken_loader(config, file):
        raise TypeError("loader invariant broken")

    monkeypatch.setattr(configuration_yaml, "load_configuration_yaml", broken_loader)

    with pytest.raises(TypeError, match="loader invariant broken"):
        configuration.load_user_preferences()


def test_user_preferences_propagate_legacy_python_error(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"
    legacy_path.write_text("raise RuntimeError('legacy config failed')\n", encoding="utf8")

    monkeypatch.setattr(configuration, "config", fresh_config())
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(
        configuration,
        "legacy_user_preferences_file",
        lambda: str(legacy_path),
    )

    with pytest.raises(RuntimeError, match="legacy config failed"):
        configuration.load_user_preferences()
