from concurrent.futures import ThreadPoolExecutor
import importlib
import sys
from threading import Barrier, Thread

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


def test_modules_import_without_materializing_values(tmp_path, monkeypatch):
    class ImportOnlyRegistry(ModuleRegistry):
        def value_from_module(self, name, module):
            raise AssertionError("modules must not materialize registry values")

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_modules",
        {"alpha": "value = object()\n"},
    )
    registry = ImportOnlyRegistry(package, "value")

    modules = registry.modules

    assert tuple(modules) == ("ALPHA",)
    assert modules["ALPHA"].__name__ == "_ase2sprkkr_registry_modules.alpha"


def test_module_import_does_not_hold_registry_lock(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_import_lock",
        {"alpha": "value = object()\n"},
    )
    registry = ModuleRegistry(package, "value")
    original_import = importlib.import_module

    def checked_import(name, package=None):
        acquired = []

        def try_registry_lock():
            locked = registry._lock.acquire(timeout=1)
            acquired.append(locked)
            if locked:
                registry._lock.release()

        thread = Thread(target=try_registry_lock)
        thread.start()
        thread.join()
        assert acquired == [True]
        return original_import(name, package)

    monkeypatch.setattr(importlib, "import_module", checked_import)

    registry["alpha"]


def test_concurrent_first_lookup_materializes_factory_once(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_concurrent",
        {
            "alpha": (
                "import time\n"
                "calls = 0\n"
                "def value():\n"
                "    global calls\n"
                "    calls += 1\n"
                "    time.sleep(0.05)\n"
                "    return object()\n"
            )
        },
    )
    registry = ModuleRegistry(package, "value", factory=True)
    workers = 8
    start = Barrier(workers)

    def lookup():
        start.wait(timeout=5)
        return registry["alpha"]

    with ThreadPoolExecutor(max_workers=workers) as executor:
        values = list(executor.map(lambda _: lookup(), range(workers)))

    assert all(value is values[0] for value in values)
    assert sys.modules["_ase2sprkkr_registry_concurrent.alpha"].calls == 1


def test_registry_lookup_errors_keep_their_categories(tmp_path, monkeypatch):
    class BrokenRegistry(ModuleRegistry):
        def prepare_value(self, name, module, value):
            raise KeyError("programming error")

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_lookup_errors",
        {
            "alpha": "value = object()\n",
            "missing_export": "other = object()\n",
            "import_failure": "import _ase2sprkkr_missing_dependency\n",
        },
    )
    registry = ModuleRegistry(package, "value")

    with pytest.raises(RegistryEntryNotFound):
        registry["missing"]

    assert "missing_export" in registry
    with pytest.raises(RegistryDefinitionError, match="required export"):
        registry["missing_export"]

    with pytest.raises(ModuleNotFoundError) as exc:
        registry["import_failure"]
    assert exc.value.name == "_ase2sprkkr_missing_dependency"

    with pytest.raises(KeyError, match="programming error"):
        BrokenRegistry(package, "value")["alpha"]


def test_registry_name_and_alias_validation(tmp_path, monkeypatch):
    class CollidingRegistry(ModuleRegistry):
        @staticmethod
        def normalize_name(name):
            return "COLLISION"

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_alias_contracts",
        {"alpha": "value = object()\n", "beta": "value = object()\n"},
    )

    with pytest.raises(RegistryDefinitionError, match="Duplicate module registry key"):
        CollidingRegistry(package, "value").names()

    with pytest.raises(RegistryDefinitionError, match="Duplicate normalized registry alias"):
        ModuleRegistry(
            package, "value",
            aliases={"first": "alpha", "FIRST": "beta"},
        )

    with pytest.raises(RegistryDefinitionError, match="shadows a module"):
        ModuleRegistry(
            package, "value", aliases={"ALPHA": "BETA"},
        ).names()

    with pytest.raises(RegistryDefinitionError, match="not a canonical module"):
        ModuleRegistry(
            package, "value",
            aliases={"FIRST": "SECOND", "SECOND": "ALPHA"},
        ).names()

    registry = ModuleRegistry(
        package, "value", aliases={"FIRST": "ALPHA"},
    )
    with pytest.raises(TypeError):
        registry.aliases["FIRST"] = "OTHER"


def test_tolerant_enumeration_skips_broken_modules_once(tmp_path, monkeypatch):
    from ..module_registry import RegistryLoadWarning

    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_tolerant",
        {
            "valid": "value = 42\n",
            "missing_export": "other = 0\n",
            "import_failure": "import _missing_registry_dependency\n",
            "factory_failure": (
                "def value():\n"
                "    raise ValueError('broken factory')\n"
            ),
        },
    )
    registry = ModuleRegistry(package, "value", factory=True)

    assert set(registry.names()) == {
        "VALID", "MISSING_EXPORT", "IMPORT_FAILURE", "FACTORY_FAILURE"
    }
    assert "_ase2sprkkr_registry_tolerant.valid" not in sys.modules
    with pytest.warns(RegistryLoadWarning) as recorded:
        assert list(registry.iter_valid_values()) == [42]
    assert len(recorded) == 3
    messages = "\n".join(str(item.message) for item in recorded)
    assert "MISSING_EXPORT" in messages.upper()
    assert "IMPORT_FAILURE" in messages.upper()
    assert "FACTORY_FAILURE" in messages.upper()

    import warnings
    with warnings.catch_warnings(record=True) as repeated:
        warnings.simplefilter("always")
        assert list(registry.iter_valid_values()) == [42]
    assert not repeated

    with pytest.raises(RegistryDefinitionError, match="required export"):
        registry["MISSING_EXPORT"]
    with pytest.raises(ModuleNotFoundError):
        registry["IMPORT_FAILURE"]
    with pytest.raises(ValueError, match="broken factory"):
        registry["FACTORY_FAILURE"]


def test_filtered_tolerant_enumeration_does_not_import_unrelated_modules(
    tmp_path, monkeypatch
):
    package = _package(
        tmp_path, monkeypatch, "_ase2sprkkr_registry_filtered",
        {
            "valid": "value = 7\n",
            "unrelated": "raise RuntimeError('unexpected import')\n",
        },
    )
    registry = ModuleRegistry(package, "value")
    assert list(registry.iter_valid_values(names=["VALID"])) == [7]
    assert "_ase2sprkkr_registry_filtered.unrelated" not in sys.modules


def test_underscore_modules_are_private_even_with_custom_include(tmp_path, monkeypatch):
    package = _package(
        tmp_path, monkeypatch, "_ase2sprkkr_registry_private",
        {"alpha": "value = 42\n", "_helper": "raise AssertionError('private module imported')\n"},
    )
    registry = ModuleRegistry(package, "value", include=lambda name: True)
    assert registry.names() == ("ALPHA",)
    assert "_HELPER" not in registry
    with pytest.raises(RegistryEntryNotFound):
        registry["_helper"]
    assert registry["alpha"] == 42
    assert "_ase2sprkkr_registry_private._helper" not in sys.modules
