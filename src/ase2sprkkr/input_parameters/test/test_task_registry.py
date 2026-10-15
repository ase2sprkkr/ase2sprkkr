import importlib

import pytest

from ...common.module_registry import RegistryDefinitionError
from ..input_parameters import InputParameters
from ..task_definition import TaskRegistry


def _external_task_package(tmp_path, monkeypatch, name, modules):
    package_root = tmp_path / f"{name}_package_root"
    package_dir = package_root / name
    package_dir.mkdir(parents=True)
    (package_dir / "__init__.py").write_text("")

    user_dir = tmp_path / f"{name}_user_tasks"
    user_dir.mkdir()
    for module_name, source in modules.items():
        (user_dir / f"{module_name}.py").write_text(source)

    monkeypatch.syspath_prepend(str(package_root))
    package = importlib.import_module(name)
    package.__path__.append(str(user_dir))
    return package


def test_task_registry_discovers_only_task_modules_without_loading_them():
    registry = InputParameters.task_registry

    names = registry.names(include_aliases=False)

    assert "SCF" in names
    assert "SECTIONS" not in names
    assert "BSFEK" not in names
    assert "BSFKK" not in names


def test_task_registry_aliases_are_task_definitions():
    canonical = InputParameters.task_registry["BSF"]
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek is not canonical
    assert bsfkk is not canonical
    assert bsfek.input_parameters_definition is canonical.input_parameters_definition
    assert bsfkk.input_parameters_definition is canonical.input_parameters_definition
    assert bsfek.name == "BSFEK"
    assert bsfkk.name == "BSFKK"


def test_task_registry_discovers_user_task_path(tmp_path, monkeypatch):
    package = _external_task_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_external_task_definitions",
        {
            "user": (
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "class _Definition:\n"
                "    name = 'USER'\n"
                "    def description(self, verbose='all'):\n"
                "        return 'User-defined task schema'\n"
                "task_definition = TaskDefinition(_Definition())\n"
            )
        },
    )
    registry = TaskRegistry(package)

    assert "USER" in registry.names(include_aliases=False)

    task = registry["USER"]

    assert task.name == "USER"
    assert task.input_parameters_definition.name == "USER"


def test_legacy_user_input_parameters_export_is_not_silently_accepted(
    tmp_path, monkeypatch
):
    package = _external_task_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_legacy_task_definitions",
        {"legacy": "input_parameters = object()\n"},
    )
    registry = TaskRegistry(package)

    assert "LEGACY" in registry
    with pytest.raises(RegistryDefinitionError, match="task_definition"):
        registry["LEGACY"]
