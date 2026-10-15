import importlib
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from ..module_registry import (
    ModuleRegistry,
    RegistryDefinitionError,
    RegistryEntryNotFound,
)


def _package(tmp_path, monkeypatch, name, modules):
    root = tmp_path / name
    root.mkdir()
    (root / "__init__.py").write_text("")
    for module_name, source in modules.items():
        (root / f"{module_name}.py").write_text(source)
    monkeypatch.syspath_prepend(str(tmp_path))
    return importlib.import_module(name)


def test_discovery_and_membership_do_not_import_modules(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_lazy",
        {
            "alpha": "value = object()\n",
            "beta": "value = object()\n",
            "helper": "raise RuntimeError('must not be imported')\n",
        },
    )
    registry = ModuleRegistry(
        package,
        "value",
        aliases={"FIRST": "ALPHA"},
        include=lambda name: name != "helper",
    )

    assert registry.names() == ("ALPHA", "BETA", "FIRST")
    assert "alpha" in registry
    assert "first" in registry
    assert "helper" not in registry
    assert "_ase2sprkkr_registry_lazy.alpha" not in sys.modules
    assert "_ase2sprkkr_registry_lazy.beta" not in sys.modules
    assert "_ase2sprkkr_registry_lazy.helper" not in sys.modules

    alpha = registry["alpha"]
    assert "_ase2sprkkr_registry_lazy.alpha" in sys.modules
    assert registry["FIRST"] is alpha
    assert "_ase2sprkkr_registry_lazy.beta" not in sys.modules


def test_base_registry_does_not_call_callable_values(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_callable",
        {"alpha": "def value():\n    return 'called'\n"},
    )
    registry = ModuleRegistry(package, "value")

    value = registry["alpha"]

    assert callable(value)
    assert value() == "called"


def test_missing_module_has_specific_lookup_error(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_missing",
        {"alpha": "value = object()\n"},
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(RegistryEntryNotFound):
        registry["missing"]


def test_malformed_module_is_not_reported_as_missing(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_missing_export",
        {"alpha": "other = object()\n"},
    )
    registry = ModuleRegistry(package, "value")

    assert "alpha" in registry
    with pytest.raises(RegistryDefinitionError, match="required export"):
        registry["alpha"]


def test_registry_does_not_hide_import_errors(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_broken",
        {"broken": "import _ase2sprkkr_missing_dependency\n"},
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(ModuleNotFoundError) as exc:
        registry["broken"]

    assert exc.value.name == "_ase2sprkkr_missing_dependency"


def test_registry_does_not_hide_value_errors(tmp_path, monkeypatch):
    class BrokenRegistry(ModuleRegistry):
        def prepare_value(self, name, module, value):
            raise KeyError("programming error")

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_programming_error",
        {"alpha": "value = object()\n"},
    )
    registry = BrokenRegistry(package, "value")

    with pytest.raises(KeyError, match="programming error"):
        registry["alpha"]


def test_duplicate_normalized_module_names_fail_explicitly(tmp_path, monkeypatch):
    class CollidingRegistry(ModuleRegistry):
        @staticmethod
        def normalize_name(name):
            return "COLLISION"

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_duplicates",
        {"alpha": "value = 1\n", "beta": "value = 2\n"},
    )
    registry = CollidingRegistry(package, "value")

    with pytest.raises(RegistryDefinitionError, match="Duplicate module registry key"):
        registry.names()


def test_duplicate_normalized_aliases_fail_explicitly(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_duplicate_aliases",
        {"alpha": "value = object()\n", "beta": "value = object()\n"},
    )

    with pytest.raises(RegistryDefinitionError, match="Duplicate normalized registry alias"):
        ModuleRegistry(
            package,
            "value",
            aliases={"first": "alpha", "FIRST": "beta"},
        )


def test_alias_may_not_shadow_a_module(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_alias_shadow",
        {"alpha": "value = object()\n", "beta": "value = object()\n"},
    )
    registry = ModuleRegistry(
        package,
        "value",
        aliases={"ALPHA": "BETA"},
    )

    with pytest.raises(RegistryDefinitionError, match="shadows a module"):
        registry.names()


def test_alias_must_target_canonical_module_directly(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_alias_chain",
        {"alpha": "value = object()\n"},
    )
    registry = ModuleRegistry(
        package,
        "value",
        aliases={"FIRST": "SECOND", "SECOND": "ALPHA"},
    )

    with pytest.raises(RegistryDefinitionError, match="not a canonical module"):
        registry.names()


def test_aliases_are_immutable(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_immutable_aliases",
        {"alpha": "value = object()\n"},
    )
    registry = ModuleRegistry(
        package,
        "value",
        aliases={"FIRST": "ALPHA"},
    )

    with pytest.raises(TypeError):
        registry.aliases["FIRST"] = "OTHER"


def test_value_materialization_is_thread_safe(tmp_path, monkeypatch):
    class CountingRegistry(ModuleRegistry):
        def __init__(self, package):
            super().__init__(package)
            self.calls = 0

        def value_from_module(self, name, module):
            self.calls += 1
            time.sleep(0.02)
            return object()

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_threads",
        {"alpha": ""},
    )
    registry = CountingRegistry(package)

    with ThreadPoolExecutor(max_workers=8) as executor:
        values = list(executor.map(lambda _: registry["alpha"], range(16)))

    assert registry.calls == 1
    assert all(value is values[0] for value in values)
