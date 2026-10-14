import importlib
import sys

import pytest

from ..module_registry import ModuleRegistry, RegistryEntryNotFound


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


def test_missing_entry_has_specific_lookup_error(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_missing",
        {"alpha": "value = object()\n"},
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(RegistryEntryNotFound):
        registry["missing"]


def test_missing_export_has_specific_lookup_error(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_missing_export",
        {"alpha": "other = object()\n"},
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(RegistryEntryNotFound):
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
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_duplicates",
        {"alpha": "value = 1\n", "ALPHA": "value = 2\n"},
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(ValueError, match="Duplicate module registry key"):
        registry.names()
