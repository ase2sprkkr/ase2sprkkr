import importlib
import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..input_parameters import InputParameters
from ..task_definition import InputParametersRegistry, TaskDefinition, TaskRegistry


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
    assert InputParameters.input_parameters_registry.definition_for_task(task) is schema
    assert "SCF" not in InputParameters.task_registry.extra_tasks.names()


def test_spec_tasks_share_schema_and_requested_task_is_default():
    schema = InputParameters.input_parameters_registry["ARPES"]
    spleed = InputParameters.create("SPLEED")

    assert schema._task_definition is None
    assert spleed._definition is schema
    assert spleed.requested_task_name == "SPLEED"
    assert spleed.TASK.TASK() == "SPLEED"
    assert not spleed.TASK.TASK.is_set()
    assert spleed.task_name == "SPLEED"

    spleed.TASK.TASK.set("ARPES")
    assert spleed.TASK.TASK.is_set()
    assert spleed.requested_task_name == "SPLEED"
    assert spleed.task_name == "ARPES"
    assert spleed.task_definition.name == "ARPES"


def test_bsf_is_schema_not_task_and_views_keep_schema_defaults():
    schema = InputParameters.input_parameters_registry["BSF"]
    ek = InputParameters.create("BSFEK")
    kk = InputParameters.create("BSFKK")

    assert schema._task_definition is None
    assert "BSF" not in InputParameters.task_registry
    assert ek._definition is schema
    assert kk._definition is schema
    assert ek.requested_task_name == "BSFEK"
    assert kk.requested_task_name == "BSFKK"
    assert ek.ENERGY.NE()[0] > 1
    assert kk.ENERGY.NE()[0] == 1
    assert not ek.ENERGY.NE.is_set()
    assert not kk.ENERGY.NE.is_set()
    assert ek.task_name == "BSFEK"
    assert kk.task_name == "BSFKK"

    ek.ENERGY.NE = [1]
    assert ek.requested_task_name == "BSFEK"
    assert ek.task_name == "BSFKK"

    with pytest.raises(RegistryEntryNotFound):
        InputParameters.create("BSF")
    assert InputParameters.is_task_name("BSF") is False


def test_input_and_task_registry_namespaces_are_distinct():
    input_names = set(InputParameters.input_parameters_registry.names())
    task_names = set(InputParameters.task_registry.names())

    assert {"ARPES", "BSF"} <= input_names
    assert "SPLEED" not in input_names
    assert "BSFEK" not in input_names
    assert "BSF" not in task_names
    assert {"ARPES", "SPLEED", "BSFEK", "BSFKK"} <= task_names


def test_input_parameters_registry_owns_input_parameter_operations():
    registry = InputParameters.input_parameters_registry

    for name in (
        "definition_for_task",
        "create_for_task",
        "read_for_task",
        "create_from_result",
        "files_for_task",
    ):
        assert callable(getattr(registry, name))
    assert not hasattr(registry, "task_definitions")


def test_task_registry_is_identity_resolution_only():
    registry = InputParameters.task_registry

    assert callable(registry.task_for_parameters)
    assert not hasattr(registry, "effective_task")
    for name in (
        "input_parameters_definition",
        "create_input_parameters",
        "read_input_parameters",
        "input_parameters_from_result",
        "input_files",
        "output_files",
        "cache_files",
    ):
        assert not hasattr(registry, name)


def test_schema_does_not_resolve_task_services():
    schema = InputParameters.input_parameters_registry["ARPES"]

    assert not hasattr(schema, "task_definition")
    assert not hasattr(schema, "initialize_for_task")
    assert not hasattr(schema, "update_from_result")
    assert callable(schema.initialize_from_result)
    assert schema.task_definition_name(InputParameters.create("SPLEED")) == "SPLEED"


def test_task_registry_resolves_task_from_values():
    parameters = InputParameters.create("ARPES")
    parameters.TASK.TASK = "SPLEED"

    task = InputParameters.task_registry.task_for_parameters(parameters)

    assert task is InputParameters.task_registry["SPLEED"]
    assert parameters.task_name == "SPLEED"
    assert parameters.requested_task_name == "ARPES"


def test_input_parameters_registry_snapshots_result_defaults():
    result = SimpleNamespace(
        last_iteration=SimpleNamespace(
            energy=SimpleNamespace(EMIN=lambda: -0.35)
        )
    )
    task = InputParameters.task_registry["JXC"]

    parameters = InputParameters.input_parameters_registry.create_from_result(
        task, result
    )

    assert parameters.task_name == "JXC"
    assert parameters.requested_task_name == "JXC"
    assert not hasattr(parameters, "from_result")
    assert parameters.task_definition.name == "JXC"
    assert parameters.ENERGY.EMIN() == -0.35
    assert parameters.ENERGY.EMIN.default_value == -0.35
    assert not parameters.ENERGY.EMIN.is_set()

    for copy_values in (False, True):
        copied = parameters.copy(copy_values=copy_values)
        assert copied.requested_task_name == "JXC"
        assert copied.ENERGY.EMIN() == -0.35
        assert copied.ENERGY.EMIN.default_value == -0.35
        assert not copied.ENERGY.EMIN.is_set()

    restored = pickle.loads(pickle.dumps(parameters))
    assert restored.requested_task_name == "JXC"
    assert restored.ENERGY.EMIN() == -0.35
    assert restored.ENERGY.EMIN.default_value == -0.35
    assert not restored.ENERGY.EMIN.is_set()

    parameters.ENERGY.EMIN = -0.35
    copied = parameters.copy()
    assert copied.ENERGY.EMIN() == -0.35
    assert copied.ENERGY.EMIN.is_set()

    parameters.ENERGY.EMIN = -0.5
    assert parameters.ENERGY.EMIN() == -0.5
    assert parameters.ENERGY.EMIN.is_set()


def test_jxc_without_result_keeps_normal_emin_default():
    parameters = InputParameters.create("JXC")

    assert not hasattr(parameters, "from_result")
    assert parameters.ENERGY.EMIN() == -0.2
    assert not parameters.ENERGY.EMIN.is_set()


def test_instance_default_can_be_cleared():
    parameters = InputParameters.create("JXC")
    parameters.ENERGY.EMIN.set_default(-0.35)

    assert parameters.ENERGY.EMIN() == -0.35
    assert not parameters.ENERGY.EMIN.is_set()

    parameters.ENERGY.EMIN.clear_default()

    assert parameters.ENERGY.EMIN() == -0.2
    assert not parameters.ENERGY.EMIN.is_set()


def test_instance_default_copy_respects_copy_values():
    parameters = InputParameters.create("JXC")
    parameters.ENERGY.GRID.set_default([1, 2])

    shallow = parameters.copy(copy_values=False)
    deep = parameters.copy(copy_values=True)

    assert shallow.ENERGY.GRID.default_value is parameters.ENERGY.GRID.default_value
    assert deep.ENERGY.GRID.default_value is not parameters.ENERGY.GRID.default_value
    np.testing.assert_equal(shallow.ENERGY.GRID(), parameters.ENERGY.GRID())
    np.testing.assert_equal(deep.ENERGY.GRID(), parameters.ENERGY.GRID())


def test_task_cannot_reuse_schema_bound_to_different_canonical_task():
    task = TaskDefinition("SCF", name="SPECIAL")

    with pytest.raises(RegistryDefinitionError, match="different canonical task metadata"):
        InputParameters.input_parameters_registry.definition_for_task(task)


def test_canonical_task_metadata_must_match_as_a_value():
    task = TaskDefinition("SCF", name="SCF", executable="different-executable")

    with pytest.raises(RegistryDefinitionError, match="different canonical task metadata"):
        InputParameters.input_parameters_registry.definition_for_task(task)


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

    assert input_registry.definition_for_task(task) is input_registry["SHARED"]


def test_old_input_parameters_export_is_rejected(tmp_path, monkeypatch):
    definitions_package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_test_legacy_ip_registry",
        {
            "legacy": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "def input_parameters():\n"
                "    return InputParametersDefinition('legacy', [])\n"
            ),
        },
    )

    registry = InputParametersRegistry(definitions_package)
    with pytest.raises(
        RegistryDefinitionError,
        match="does not define required export 'input_parameters_definition'",
    ):
        registry["LEGACY"]


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
