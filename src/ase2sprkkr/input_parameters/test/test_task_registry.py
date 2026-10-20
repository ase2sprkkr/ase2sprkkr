import importlib
import io
import os
import pickle
import subprocess
import sys
from collections.abc import Mapping
from types import SimpleNamespace

import numpy as np
import pytest

from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..input_parameters import InputParameters, _LazyInputParametersDefinitions
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
    task_registry = InputParameters.task_registry
    registry = InputParameters.input_parameters_registry
    arpes_task = task_registry["ARPES"]
    spleed_task = task_registry["SPLEED"]
    assert registry.definition_for_task(arpes_task) is schema
    assert registry.definition_for_task(spleed_task) is schema
    assert arpes_task.input_parameters_definition_name == schema.name.upper()
    assert spleed_task.input_parameters_definition_name == schema.name.upper()
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

def test_jxc_and_dmi_are_distinct_tasks_sharing_one_schema():
    schema = InputParameters.input_parameters_registry["JXC"]
    jxc = InputParameters.task_registry["JXC"]
    dmi = InputParameters.task_registry["DMI"]
    assert schema._task_definition is None
    assert jxc.name == "JXC"
    assert dmi.name == "DMI"
    assert dmi is not jxc
    assert jxc.input_parameters_definition_name == "JXC"
    assert dmi.input_parameters_definition_name == "JXC"
    assert InputParameters.input_parameters_registry.definition_for_task(jxc) is schema
    assert InputParameters.input_parameters_registry.definition_for_task(dmi) is schema
    assert "JXC" in InputParameters.task_registry.names()
    assert "DMI" in InputParameters.task_registry.names()
    assert InputParameters.task_registry["JXC"] is jxc
    assert InputParameters.task_registry["DMI"] is dmi


def test_shared_schema_task_lookup_is_lazy(tmp_path, monkeypatch):
    definitions = _package(
        tmp_path, monkeypatch, "_lazy_distinct_task_ipds", {
            "shared": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "def input_parameters_definition():\n"
                "    return InputParametersDefinition('shared', [])\n"
            ),
            "broken": (
                "def input_parameters_definition():\n"
                "    raise AssertionError('Unrelated IPD was constructed')\n"
            ),
        },
    )
    tasks = _package(tmp_path, monkeypatch, "_lazy_distinct_task_entries", {
        "first": (
            "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
            "def task_definition():\n"
            "    return TaskDefinition('SHARED', name='FIRST')\n"
        ),
        "second": (
            "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
            "def task_definition():\n"
            "    return TaskDefinition('SHARED', name='SECOND')\n"
        ),
    })
    ip_registry = InputParametersRegistry(definitions)
    registry = TaskRegistry(ip_registry, tasks)
    assert registry["FIRST"].name == "FIRST"
    assert "SHARED" in ip_registry._values
    assert "BROKEN" not in ip_registry._values
    assert "SECOND" not in registry.extra_tasks._values
    assert "_lazy_distinct_task_ipds.broken" not in __import__("sys").modules


def test_read_for_distinct_task_requires_matching_native_identity():
    from ...input_parameters.task_definition import TaskMismatchError
    ip_registry = InputParameters.input_parameters_registry
    dmi = InputParameters.task_registry["DMI"]
    jxc = InputParameters.task_registry["JXC"]
    params = InputParameters.create("DMI")
    params.CONTROL.POTFIL = "Fe.pot"
    text = params.to_string(validate=True)
    parsed = ip_registry.read_for_task(dmi, io.StringIO(text))
    assert parsed.task_name == "DMI"
    assert parsed.task_definition is dmi
    assert parsed.requested_task_name == "DMI"
    with pytest.raises(TaskMismatchError):
        ip_registry.read_for_task(jxc, io.StringIO(text))


def test_create_from_dmi_result_keeps_dmi_defaults():
    result = SimpleNamespace(
        last_iteration=SimpleNamespace(
            energy=SimpleNamespace(EMIN=lambda: -0.35)
        )
    )
    dmi = InputParameters.task_registry["DMI"]
    params = InputParameters.input_parameters_registry.create_from_result(dmi, result)
    assert params.task_name == "DMI"
    assert params.requested_task_name == "DMI"
    assert params.TASK.DMI() and not params.TASK.JXC()
    assert params.ENERGY.EMIN() == -0.35


def test_compatibility_definition_mapping_loads_only_requested_schema(
    tmp_path, monkeypatch
):
    package = _package(
        tmp_path,
        monkeypatch,
        "_lazy_user_definition_entries",
        {
            "valid": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "def input_parameters_definition():\n"
                "    return InputParametersDefinition('valid', [])\n"
            ),
            "marpes": "input_parameters = object()\n",
        },
    )
    mapping = _LazyInputParametersDefinitions(InputParametersRegistry(package))

    assert isinstance(mapping, Mapping)
    assert set(mapping) == {"VALID", "MARPES"}
    assert "MARPES" in mapping  # Discovery is not module materialization.
    assert mapping["VALID"].name == "valid"

    with pytest.raises(
        RegistryDefinitionError,
        match="does not define required export 'input_parameters_definition'",
    ):
        mapping["MARPES"]


def test_legacy_user_hook_can_read_arpes_without_loading_broken_user_schema(
    tmp_path,
):
    # Reproduce an import-time preference using InputParameters.definitions.
    # The separate user-installed MARPES module is intentionally malformed.
    config_dir = tmp_path / "ase2sprkkr"
    custom_schemas = config_dir / "input_parameters"
    custom_schemas.mkdir(parents=True)
    (custom_schemas / "marpes.py").write_text(
        "input_parameters = object()\n"
    )
    (config_dir / "__init__.py").write_text(
        "from ase2sprkkr.input_parameters.input_parameters import InputParameters\n"
        "InputParameters.definitions['ARPES'].mpi = False\n"
        "assert InputParameters.definitions['ARPES'].name == 'arpes'\n"
        "assert InputParameters.definitions['ARPES'].mpi is False\n"
    )
    env = os.environ.copy()
    env["XDG_CONFIG_HOME"] = str(tmp_path)
    env.pop("ASE2SPRKKR_NO_USER_PROFILE", None)

    # Do not execute in XDG_CONFIG_HOME: its ase2sprkkr directory would
    # shadow the installed package on sys.path in a fresh Python process.
    working_dir = tmp_path / "work"
    working_dir.mkdir()
    completed = subprocess.run(
        [sys.executable, "-c", "import ase2sprkkr"],
        env=env,
        cwd=working_dir,
        capture_output=True,
        text=True,
        timeout=45,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_task_enumeration_ignores_broken_schema_and_extra_task(
    tmp_path, monkeypatch
):
    from ...common.module_registry import RegistryLoadWarning

    definitions = _package(
        tmp_path, monkeypatch, "_tolerant_task_schemas", {
            "good": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "def input_parameters_definition():\n"
                "    definition = InputParametersDefinition('good', [])\n"
                "    definition._task_definition = TaskDefinition('GOOD')\n"
                "    return definition\n"
            ),
            "marpes": "input_parameters = object()\n",
            "shared": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "def input_parameters_definition():\n"
                "    return InputParametersDefinition('shared', [])\n"
            ),
        },
    )
    tasks = _package(
        tmp_path, monkeypatch, "_tolerant_extra_tasks", {
            "extra": (
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "def task_definition():\n"
                "    return TaskDefinition('SHARED', name='EXTRA')\n"
            ),
            "invalid": "task_definition = object()\n",
        },
    )
    registry = TaskRegistry(InputParametersRegistry(definitions), tasks)
    with pytest.warns(RegistryLoadWarning) as record:
        assert set(registry.names()) == {"GOOD", "EXTRA"}
    assert len(record) == 2
    assert {registry[name].name for name in ("GOOD", "EXTRA")} == {"GOOD", "EXTRA"}
    with pytest.raises(RegistryDefinitionError, match="required export"):
        registry.input_parameters_registry["MARPES"]
    with pytest.raises(RegistryDefinitionError, match="must be TaskDefinition"):
        registry.extra_tasks["INVALID"]


def test_cli_and_input_autodetection_with_broken_user_schema(tmp_path):
    config_dir = tmp_path / "ase2sprkkr"
    custom = config_dir / "input_parameters"
    custom.mkdir(parents=True)
    (custom / "marpes.py").write_text("input_parameters = object()\n")
    work = tmp_path / "work"
    work.mkdir()
    env = os.environ.copy()
    env["XDG_CONFIG_HOME"] = str(tmp_path)
    env.pop("ASE2SPRKKR_NO_USER_PROFILE", None)

    examples = os.path.join(os.path.dirname(__file__), "..", "examples", "scf.in")

    commands = (
        ([sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "--help"], False),
        ([sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "run", "--help"], False),
        ([sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "info", "--task=scf"], False),
        ([sys.executable, "-c", (
            "from ase2sprkkr.input_parameters.input_parameters import InputParameters; "
            "assert 'SCF' in InputParameters.task_registry.names(); "
            "assert InputParameters.from_file(" + repr(examples) + ").task_name == 'SCF'"
        )], True),
    )
    for command, should_warn in commands:
        completed = subprocess.run(
            command, cwd=work, env=env, capture_output=True, text=True,
            timeout=45, check=False,
        )
        assert completed.returncode == 0, (
            f"{command!r} failed: {completed.stderr}\n{completed.stdout}"
        )
        assert ("Skipping invalid registry entry" in completed.stderr) == should_warn
