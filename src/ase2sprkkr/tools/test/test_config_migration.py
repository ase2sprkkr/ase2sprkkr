from types import SimpleNamespace

from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.tools.commands import config as config_command


def test_migrate_cleans_migrated_assignments_but_keeps_python_hook(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    init_path = tmp_path / "__init__.py"
    yaml_path = tmp_path / "config.yaml"
    init_path.write_text(
        "from ase2sprkkr.configuration import config\n"
        "config.executables.dir = '/opt/sprkkr'\n"
        "config.running.mpi = 4\n"
        "custom_python_hook = 42\n",
        encoding="utf8",
    )
    args = SimpleNamespace(
        keep_init=False,
        clean_init=False,
        delete_init=False,
    )

    config_command._migrate(args, str(init_path), str(yaml_path))

    assert yaml_path.is_file()
    assert init_path.is_file()
    text = init_path.read_text(encoding="utf8")
    assert "config.executables.dir" not in text
    assert "config.running.mpi" not in text
    assert "from ase2sprkkr.configuration import config" not in text
    assert "custom_python_hook = 42" in text
