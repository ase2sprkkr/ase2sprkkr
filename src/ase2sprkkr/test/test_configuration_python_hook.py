import ase2sprkkr.configuration as configuration


def test_python_configuration_hook_is_loaded_with_yaml(tmp_path, monkeypatch):
    yaml_path = tmp_path / "config.yaml"
    init_path = tmp_path / "__init__.py"
    yaml_path.write_text(
        "default_profile: yaml\nprofiles:\n  yaml: {}\n",
        encoding="utf8",
    )
    init_path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.profiles.add('from-python')\n"
        "config.default_profile = 'from-python'\n",
        encoding="utf8",
    )

    config = configuration.definition.create_object()
    monkeypatch.setattr(configuration, "config", config)
    monkeypatch.setattr(configuration, "user_preferences_file", lambda: str(yaml_path))
    monkeypatch.setattr(
        configuration, "legacy_user_preferences_file", lambda: str(init_path)
    )

    configuration.load_user_preferences()

    assert "from-python" in config.profiles
    assert "yaml" in config.profiles
    assert config.default_profile() == "yaml"
