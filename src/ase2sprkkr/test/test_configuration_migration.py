import os

import pytest

from ase2sprkkr import configuration
from ase2sprkkr.configuration import definition
from ase2sprkkr.configuration_migration import (
    clean_legacy_configuration,
    load_legacy_configuration,
    migrate_legacy_configuration,
)


@pytest.mark.parametrize("config_import", [
    "from ase2sprkkr import config",
    "from ase2sprkkr.configuration import config",
    "import ase2sprkkr; config = ase2sprkkr.config",
])
def test_migration_isolates_all_public_config_imports(tmp_path, config_import):
    before = configuration.config.running.print_output()
    requested = before is False
    source = tmp_path / "__init__.py"
    source.write_text(
        f"{config_import}\n"
        f"import os; assert os.getpid() != {os.getpid()}\n"
        f"config.running.print_output = {requested!r}\n"
    )

    migrated = load_legacy_configuration(source, definition)
    assert migrated.running.print_output() is requested
    assert configuration.config.running.print_output() == before


def test_failed_migration_cannot_mutate_parent_configuration(tmp_path):
    before = configuration.config.running.print_output()
    requested = before is False
    source = tmp_path / "__init__.py"
    source.write_text(
        "from ase2sprkkr import config\n"
        f"config.running.print_output = {requested!r}\n"
        "raise ValueError('injected hook failure')\n"
    )
    with pytest.raises(RuntimeError, match="injected hook failure"):
        load_legacy_configuration(source, definition)
    assert configuration.config.running.print_output() == before


def test_clean_legacy_configuration_removes_comments_and_extra_blank_lines(tmp_path):
    path = tmp_path / "__init__.py"
    path.write_text(
        "# ASE2SRPKKR configuration file\n"
        "# -------------------------------\n"
        "\n"
        "# Please, DO NOT comment the following line.\n"
        "from ase2sprkkr.configuration import config\n"
        "from ase2sprkkr import InputParameters\n"
        "InputParameters.definitions['ARPES'].mpi = False\n"
        "\n"
        "# This file is pure python and it is executed when ase2sprkkr is imported.\n"
        "# Place it into the following path:\n"
        "# /home/user/.config/ase2sprkkr/__init__.py\n"
        "\n"
        "# This string is appended to the runned executables\n"
        'config.executables.suffix = "" #"8.6"\n'
        "\n"
        "# Do you want to run the executables from a specific directory?\n"
        'config.executables.dir = "/opt/sprkkr"\n'
        "\n"
        "# Uncomment, if you don't want to run empty-spheres finding by default\n"
        "config.running.empty_spheres = False\n"
        "\n"
        "# Set to False if MPI should not be used. Or set to the number of processor,\n"
        "# or just command line to run mpi programm, e.g.: [ '/usr/bin/mpiexec', '-n', '4' ]\n"
        "# config.running.mpi = []\n"
        "\n"
        "# You can change the verbosity of the output by setting the following to False or True\n"
        "# config.running.print_output = 'info'\n"
        "\n"
        "# Authentication token to Nomad.\n"
        'config.nomad.token = "TOKEN"\n',
        encoding="utf8",
    )

    assert clean_legacy_configuration(path) == 4
    assert path.read_text(encoding="utf8") == (
        "# ASE2SRPKKR configuration file\n"
        "# -------------------------------\n"
        "\n"
        "# Please, DO NOT comment the following line.\n"
        "from ase2sprkkr import InputParameters\n"
        "InputParameters.definitions['ARPES'].mpi = False\n"
        "\n"
        "# This file is pure python and it is executed when ase2sprkkr is imported.\n"
        "# Place it into the following path:\n"
        "# /home/user/.config/ase2sprkkr/__init__.py\n"
        "\n"
        "# Set to False if MPI should not be used. Or set to the number of processor,\n"
        "# or just command line to run mpi programm, e.g.: [ '/usr/bin/mpiexec', '-n', '4' ]\n"
        "# config.running.mpi = []\n"
        "\n"
        "# You can change the verbosity of the output by setting the following to False or True\n"
        "# config.running.print_output = 'info'\n"
        "\n"
    )


def test_migrate_legacy_configuration_backs_up_original_file(tmp_path):
    legacy_path = tmp_path / "__init__.py"
    yaml_path = tmp_path / "config.yaml"
    original = (
        "from ase2sprkkr.configuration import config\n"
        "config.executables.dir = '/opt/legacy'\n"
        "config.running.mpi = 3\n"
    )
    legacy_path.write_text(original, encoding="utf8")

    migrate_legacy_configuration(legacy_path, yaml_path, definition)

    backup = tmp_path / "__init__.py.bak"
    assert backup.read_text(encoding="utf8") == original
    assert legacy_path.read_text(encoding="utf8") == original
    assert yaml_path.is_file()
