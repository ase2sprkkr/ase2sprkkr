import pytest

import ase2sprkkr.configuration as configuration
import ase2sprkkr.configuration_yaml as configuration_yaml
from ase2sprkkr.configuration import definition


def fresh_config():
    return definition.create_object()


@pytest.mark.parametrize("source", ("yaml", "legacy-python"))
def test_user_preferences_propagate_loader_errors(tmp_path, monkeypatch, source):
    yaml_path = tmp_path / "config.yaml"
    legacy_path = tmp_path / "__init__.py"

    monkeypatch.setattr(configuration, "config", fresh_config())
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(
        configuration, "legacy_user_preferences_file", lambda: str(legacy_path),
    )

    if source == "yaml":
        yaml_path.write_text("profiles: {}\n", encoding="utf8")

        def broken_loader(config, file):
            raise TypeError("loader invariant broken")

        monkeypatch.setattr(configuration_yaml, "load_configuration_yaml", broken_loader)
        exception, message = TypeError, "loader invariant broken"
    else:
        legacy_path.write_text("raise RuntimeError('legacy config failed')\n", encoding="utf8")
        exception, message = RuntimeError, "legacy config failed"

    with pytest.raises(exception, match=message):
        configuration.load_user_preferences()


@pytest.mark.parametrize("bad_yaml", (
    "profiles: [\n", "- not a mapping\n",
))
def test_unparseable_user_yaml_warns_but_preserves_file(
    tmp_path, monkeypatch, bad_yaml,
):
    from ase2sprkkr.configuration_yaml import (
        YamlConfigurationWarning, InvalidYamlDocumentError,
    )
    from ruamel.yaml.error import YAMLError

    path = tmp_path / "config.yaml"
    path.write_text(bad_yaml, encoding="utf8")
    monkeypatch.setattr(configuration, "config", fresh_config())
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(path))
    monkeypatch.setattr(
        configuration, "legacy_user_preferences_file",
        lambda: str(tmp_path / "__init__.py"),
    )

    with pytest.warns(YamlConfigurationWarning, match="were not loaded"):
        configuration.load_user_preferences()
    assert configuration.config.default_profile() is None
    assert path.read_text(encoding="utf8") == bad_yaml

    # Saving must not overwrite an invalid document with default values.
    with pytest.raises((YAMLError, InvalidYamlDocumentError)):
        configuration.config.save_user_preferences()
    assert path.read_text(encoding="utf8") == bad_yaml

    # Correcting the file clears the reported error and restores normal load.
    path.write_text("running:\n  empty_spheres: true\n", encoding="utf8")
    configuration.load_user_preferences()
    assert configuration.config.running.empty_spheres() is True
