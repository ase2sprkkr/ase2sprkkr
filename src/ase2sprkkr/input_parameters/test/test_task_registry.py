import importlib
import io
import os
import subprocess
import sys
import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from ...common.module_registry import RegistryDefinitionError, RegistryEntryNotFound
from ..input_parameters import InputParameters
from ..task_definition import InputParametersRegistry, TaskDefinition, TaskRegistry
from ...common.warnings import DataValidityWarning


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


def test_spec_tasks_have_specialized_schemas_and_requested_task_is_default():
    schema = InputParameters.input_parameters_registry["ARPES"]
    spleed = InputParameters.create("SPLEED")

    assert schema._task_definition is None
    assert spleed._definition is InputParameters.input_parameters_registry["SPLEED"]
    assert spleed._definition is not schema
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

    with pytest.warns(DataValidityWarning, match='ENERGY.NE selects BSF-KK, but BSF-EK was requested'):
        ek.ENERGY.NE = [1]
    assert ek.requested_task_name == "BSFEK"
    assert ek.task_name == "BSFKK"

    assert InputParameters.for_definition("BSF")._definition is schema
    assert InputParameters.create("BSF")._definition is schema
    with pytest.raises(RegistryEntryNotFound):
        InputParameters.for_task("BSF")
    assert InputParameters.is_task_name("BSF") is False


def test_input_and_task_registry_namespaces_are_distinct():
    input_names = set(InputParameters.input_parameters_registry.names())
    task_names = set(InputParameters.task_registry.names())

    assert {"ARPES", "BSF"} <= input_names
    assert {"SPLEED", "AIPES", "BAND"} <= input_names
    assert "BSFEK" not in input_names
    assert "BSF" not in task_names
    assert {"ARPES", "SPLEED", "BSFEK", "BSFKK"} <= task_names


def test_input_parameters_factory_entry_points():
    scf_schema = InputParameters.definition("SCF")

    from_definition = InputParameters.for_definition("SCF")
    assert from_definition._definition is scf_schema
    assert from_definition.requested_task_name is None

    from_task = InputParameters.for_task("SCF")
    assert from_task._definition is scf_schema
    assert from_task.requested_task_name == "SCF"
    assert InputParameters.create("SCF").requested_task_name == "SCF"

    # The schema and task namespaces are deliberately different.
    bsf_schema = InputParameters.definition("BSF")
    assert InputParameters.for_definition("BSF")._definition is bsf_schema
    assert InputParameters.create("BSF")._definition is bsf_schema
    with pytest.raises(RegistryEntryNotFound):
        InputParameters.for_task("BSF")
    with pytest.raises(RegistryEntryNotFound):
        InputParameters.for_definition("BSFEK")

    parameters = InputParameters.create("BSFEK")
    assert parameters.requested_task_name == "BSFEK"
    assert InputParameters.create(parameters) is parameters

    filename = os.path.join(os.path.dirname(__file__), "..", "examples", "scf.in")
    assert InputParameters.create(filename).task_name == "SCF"
    with open(filename) as stream:
        assert InputParameters.create(stream).task_name == "SCF"


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


def test_explicit_task_can_reuse_schema_with_noncanonical_name():
    task = TaskDefinition("SCF", name="SPECIAL")
    registry = InputParameters.input_parameters_registry
    assert registry.definition_for_task(task) is registry["SCF"]


def test_explicit_task_can_override_canonical_executable():
    task = TaskDefinition("SCF", name="SCF", executable="different-executable")
    registry = InputParameters.input_parameters_registry
    assert registry.definition_for_task(task) is registry["SCF"]


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




def test_expected_task_does_not_supply_missing_shared_task_selector():
    """An expected task only fills a missing selector in incomplete mode."""
    from pathlib import Path
    import pyparsing as pp
    from ...common.warnings import DataValidityWarning
    from ..task_definition import TaskMismatchError

    example = Path(__file__).resolve().parent.parent / "examples" / "arpes.in"
    original = example.read_text()
    incomplete = original.replace("TASK     ARPES\n", "TASK\n", 1)
    assert incomplete != original

    registry = InputParameters.input_parameters_registry
    spleed = InputParameters.task_registry["SPLEED"]

    # Strict parsing still requires the native task selector.
    with pytest.raises(pp.ParseBaseException, match="TASK"):
        registry.read_for_task(spleed, io.StringIO(incomplete))

    # Incomplete parsing is an explicit opt-in, with a diagnostic.
    with pytest.warns(DataValidityWarning, match="TASK.TASK"):
        parsed = registry.read_for_task(
            spleed, io.StringIO(incomplete), allow_incomplete=True
        )
    assert parsed.task_name == "SPLEED"
    assert parsed.requested_task_name == "SPLEED"
    assert not parsed.TASK.TASK.is_set()

    # An explicitly supplied, incompatible selector cannot be overridden.
    with pytest.raises(TaskMismatchError):
        registry.read_for_task(
            spleed, io.StringIO(original), allow_incomplete=True
        )


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


def _optional_task_file(parameters):
    return "{dataset}.optional" if parameters.CONTROL.PRINT() else None


def test_file_template_callbacks_skip_none_and_support_pickle(tmp_path):
    parameters = InputParameters.create("ARPES")
    parameters.CONTROL.DATASET = "Fe"
    parameters.CONTROL.PRINT = 0
    task = TaskDefinition(
        "ARPES", name="FILE_CALLBACK_TEST",
        input_files=("{dataset}.required", _optional_task_file),
    )
    registry = InputParameters.input_parameters_registry
    assert registry.files_for_task(task, parameters, str(tmp_path), "input") == (
        str(tmp_path / "Fe.required"),
    )
    parameters.CONTROL.PRINT = 1
    assert registry.files_for_task(task, parameters, str(tmp_path), "input") == (
        str(tmp_path / "Fe.required"), str(tmp_path / "Fe.optional"),
    )
    restored = pickle.loads(pickle.dumps(task))
    assert restored == task
    assert restored.input_files[1].pathname is _optional_task_file




def test_definitions_mapping_does_not_load_unrelated_definitions(monkeypatch):
    from ase2sprkkr.input_parameters.input_parameters import (
        _LazyInputParametersDefinitions,
    )

    class Registry:
        def __init__(self):
            self.loaded = []

        def names(self):
            return ("SCF", "BROKEN")

        def __getitem__(self, name):
            self.loaded.append(name)
            if name == "BROKEN":
                raise ImportError("broken user-provided definition")
            return "scf-schema"

    registry = Registry()
    mapping = _LazyInputParametersDefinitions(registry)
    assert list(mapping) == ["SCF", "BROKEN"]
    assert len(mapping) == 2
    assert registry.loaded == []
    assert mapping["SCF"] == "scf-schema"
    assert registry.loaded == ["SCF"]
    with pytest.raises(ImportError, match="broken user-provided"):
        mapping["BROKEN"]


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
    assert {task.name for task in registry.values()} == {"GOOD", "EXTRA"}
    assert {registry[name].name for name in ("GOOD", "EXTRA")} == {"GOOD", "EXTRA"}
    with pytest.raises(RegistryDefinitionError, match="required export"):
        registry.input_parameters_registry["MARPES"]
    with pytest.raises(RegistryDefinitionError, match="must be TaskDefinition"):
        registry.extra_tasks["INVALID"]


def test_task_enumeration_detects_valid_name_collision(tmp_path, monkeypatch):
    definitions = _package(
        tmp_path, monkeypatch, "_duplicate_task_schema", {
            "same": (
                "from ase2sprkkr.input_parameters.input_parameters_definitions "
                "import InputParametersDefinition\n"
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "def input_parameters_definition():\n"
                "    definition = InputParametersDefinition('same', [])\n"
                "    definition._task_definition = TaskDefinition('SAME')\n"
                "    return definition\n"
            ),
        },
    )
    tasks = _package(
        tmp_path, monkeypatch, "_duplicate_task_extras", {
            "same": (
                "from ase2sprkkr.input_parameters.task_definition import TaskDefinition\n"
                "def task_definition():\n"
                "    return TaskDefinition('SAME', name='SAME')\n"
            ),
        },
    )
    registry = TaskRegistry(InputParametersRegistry(definitions), tasks)
    with pytest.raises(RegistryDefinitionError, match="defined both"):
        registry.names()
    with pytest.raises(RegistryDefinitionError, match="defined both"):
        registry["SAME"]


def test_user_preset_ipd_reuses_canonical_task_without_autodetection_ambiguity(tmp_path):
    """An IPD copy has its own registry key, but does not declare a new task."""
    config_dir = tmp_path / "ase2sprkkr" / "input_parameters"
    config_dir.mkdir(parents=True)
    (config_dir / "myarpes.py").write_text(
        "from ase2sprkkr.input_parameters.input_parameters import InputParameters\n"
        "def input_parameters_definition():\n"
        "    definition = InputParameters.definition('ARPES').copy()\n"
        "    definition['TAU']['NKTAB'].default_value = 555\n"
        "    return definition\n"
    )
    (config_dir / "myscf.py").write_text(
        "from ase2sprkkr.input_parameters.input_parameters import InputParameters\n"
        "def input_parameters_definition():\n"
        "    return InputParameters.definition('SCF').copy()\n"
    )
    env = os.environ.copy()
    env["XDG_CONFIG_HOME"] = str(tmp_path)
    env.pop("ASE2SPRKKR_NO_USER_PROFILE", None)
    arpes_example = os.path.join(os.path.dirname(__file__), "..", "examples", "arpes.in")
    source = (
        "from ase2sprkkr.input_parameters.input_parameters import InputParameters as IP\n"
        "arpes = IP.for_definition('MYARPES')\n"
        "assert IP.create('MYARPES').TAU.NKTAB() == 555\n"
        "assert arpes._definition is IP.definition('MYARPES')\n"
        "assert arpes.task_name == 'ARPES'\n"
        "assert arpes.TAU.NKTAB() == 555\n"
        "assert arpes.task_definition_required() is IP.task_registry['ARPES']\n"
        "scf = IP.for_definition('MYSCF')\n"
        "assert scf.task_definition_required() is IP.task_registry['SCF']\n"
        "assert 'MYARPES' not in IP.task_registry.names()\n"
        "assert 'MYSCF' not in IP.task_registry.names()\n"
        f"file = {arpes_example!r}\n"
        "assert IP.from_file(file)._definition is IP.definition('ARPES')\n"
        "parsed = IP.input_parameters_registry.read_for_task(\n"
        "    IP.task_registry['ARPES'], file, input_definition='MYARPES')\n"
        "assert parsed._definition is IP.definition('MYARPES')\n"
        "assert parsed.task_name == 'ARPES'\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", source],
        env=env, capture_output=True, text=True, timeout=45, check=False,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout

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



@pytest.mark.parametrize("name", ("ARPES", "AIPES", "SPLEED", "BAND"))
def test_spec_input_requires_explicit_task_selector(name):
    import pyparsing as pp

    task = InputParameters.task_registry[name]
    definition = InputParameters.definition(name)
    parameters = InputParameters.create(name)
    parameters.CONTROL.POTFIL = "Fe.pot"
    text = parameters.to_string(validate=True)
    assert "TASK" in text

    # A present TASK section without its positional selector must be rejected.
    missing_selector = text.replace("TASK     " + name, "TASK     ", 1)
    if missing_selector == text:
        import re
        missing_selector = re.sub(
            rf"(?m)^(TASK\s+){name}\b", r"\1", text, count=1
        )
    assert missing_selector != text
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string(missing_selector)
    with pytest.raises(pp.ParseBaseException):
        InputParameters.input_parameters_registry.read_for_task(
            task, io.StringIO(missing_selector)
        )
