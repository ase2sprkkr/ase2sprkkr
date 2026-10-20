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
