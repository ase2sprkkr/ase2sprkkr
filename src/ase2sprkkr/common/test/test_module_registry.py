import importlib
import sys

import pytest

from ..module_registry import ModuleRegistry


def _package(tmp_path, monkeypatch, name, modules):
    root = tmp_path / name
    root.mkdir()
    (root / "__init__.py").write_text("")
    for module_name, source in modules.items():
        (root / f"{module_name}.py").write_text(source)
    monkeypatch.syspath_prepend(str(tmp_path))
    return importlib.import_module(name)


def test_registry_discovers_without_importing_and_loads_lazily(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_lazy",
        {
            "alpha": "value = object()\n",
            "beta": "value = object()\n",
            "helper": "raise RuntimeError('helper must not be imported')\n",
        },
    )
    registry = ModuleRegistry(
        package,
        lambda module: module.value,
        normalize=str.upper,
        include=lambda name: name != "helper",
        aliases={"FIRST": "ALPHA"},
    )

    assert tuple(registry) == ("ALPHA", "BETA")
    assert "_ase2sprkkr_registry_lazy.alpha" not in sys.modules
    assert "alpha" in registry
    assert "first" in registry
    assert "helper" not in registry

    alpha = registry["alpha"]
    assert "_ase2sprkkr_registry_lazy.alpha" in sys.modules
    assert registry["ALPHA"] is alpha
    assert registry["FIRST"] is alpha
    assert "_ase2sprkkr_registry_lazy.beta" not in sys.modules


def test_registry_does_not_hide_import_errors(tmp_path, monkeypatch):
    package = _package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_registry_broken",
        {"broken": "import _ase2sprkkr_missing_dependency\n"},
    )
    registry = ModuleRegistry(package, lambda module: module)

    with pytest.raises(ModuleNotFoundError) as exc:
        registry["broken"]

    assert exc.value.name == "_ase2sprkkr_missing_dependency"


def test_registry_requires_a_package():
    import math

    with pytest.raises(TypeError, match="is not a package"):
        ModuleRegistry(math, lambda module: module)
