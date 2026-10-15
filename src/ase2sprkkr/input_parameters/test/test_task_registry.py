import importlib

import pytest

from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..input_parameters import InputParameters
from ..task_definition import InputParametersRegistry, TaskRegistry


def _package(tmp_path, monkeypatch, name, modules):
    root = tmp_path / name
    root.mkdir()
    (root / "__init__.py").write_text("")
    for module_name, source in modules.items():
        (root / f"{module_name}.py").write_text(source)
    monkeypatch.syspath_prepend(str(tmp_path))
    return importlib.import_module(name)


def test_standard_task_comes_from_input_parameters_definition():
    schema = InputParameters.input_parameters_registry["SCF"]
    task = InputParameters.task_registry["SCF"]

    assert schema._task_definition is task
    assert task.input_parameters_definition_name == "SCF"
    assert task.input_parameters_definition is schema
    assert "SCF" not in InputParameters.task_registry.extra_tasks.names()


def test_spec_tasks_share_schema_and_requested_task_is_materialized():
    schema = InputParameters.input_parameters_registry["ARPES"]
    spleed = InputParameters.create("SPLEED")

    assert schema._task_definition is None
    assert spleed._definition is schema
    assert spleed.TASK.TASK() == "SPLEED"
    assert spleed.task_definition.name == "SPLEED"

    spleed.TASK.TASK.set("ARPES")
    assert spleed.task_name == "SPLEED"
    assert spleed.task_definition.name == "ARPES"


def test_bsf_is_schema_not_task_and_views_keep_schema_defaults():
    schema = InputParameters.input_parameters_registry["BSF"]
    ek = InputParameters.create("BSFEK")
    kk = InputParameters.create("BSFKK")

    assert schema._task_definition is None
    assert "BSF" not in InputParameters.task_registry
    assert ek._definition is schema
    assert kk._definition is schema
    assert ek.ENERGY.NE()[0] > 1
    assert kk.ENERGY.NE()[0] == 1
    assert ek.task_definition.name == "BSFEK"
    assert kk.task_definition.name == "BSFKK"

    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create("BSF")

    generic = InputParameters.create_input_parameters("BSF")
    assert generic._definition is schema


def test_input_and_task_registry_namespaces_are_distinct():
    input_names = set(InputParameters.input_parameters_registry.names())
    task_names = set(InputParameters.task_registry.names())

    assert {"ARPES", "BSF"} <= input_names
    assert "SPLEED" not in input_names
    assert "BSFEK" not in input_names
    assert "BSF" not in task_names
    assert {"ARPES", "SPLEED", "BSFEK", "BSFKK"} <= task_names


def test_extra_task_resolves_schema_from_own_input_registry(tmp_path, monkeypatch):
    definitions_package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_test_ip_registry",
        {
            "shared": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "def input_parameters_definition():\n"
                "    return InputParametersDefinition('shared', [])\n"
            ),
        },
    )
    tasks_package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_test_task_registry",
        {
            "special": (
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "def task_definition():\n"
                "    return TaskDefinition('SHARED', name='SPECIAL')\n"
            ),
        },
    )

    input_registry = InputParametersRegistry(definitions_package)
    task_registry = TaskRegistry(input_registry, tasks_package)
    task = task_registry["SPECIAL"]

    assert task.input_parameters_definition is input_registry["SHARED"]


def test_malformed_extra_task_is_definition_error(tmp_path, monkeypatch):
    definitions_package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_test_empty_ip_registry",
        {},
    )
    tasks_package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_test_bad_task_registry",
        {"broken": "task_definition = object()\n"},
    )

    registry = TaskRegistry(InputParametersRegistry(definitions_package), tasks_package)
    with pytest.raises(RegistryDefinitionError, match="must be TaskDefinition"):
        registry["BROKEN"]
