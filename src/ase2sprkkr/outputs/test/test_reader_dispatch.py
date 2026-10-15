import importlib
import sys
from types import SimpleNamespace

import pytest

from ...common.module_registry import RegistryDefinitionError
from ...input_parameters.task_definition import TaskDefinition
from .. import reader_registry as reader_registry_module
from ..reader_registry import OutputReaderRegistry, reader_registry
from ..readers.default import DefaultOutputReader
from ..task_result import KkrOutputReader


def _reader_package(tmp_path, monkeypatch, name, modules):
    root = tmp_path / name
    root.mkdir()
    (root / "__init__.py").write_text("")
    for module_name, source in modules.items():
        (root / f"{module_name}.py").write_text(source)
    monkeypatch.syspath_prepend(str(tmp_path))
    return importlib.import_module(name)


def test_missing_task_reader_uses_default_reader():
    assert reader_registry.class_for_task("definitely_missing_reader") is DefaultOutputReader


def test_reader_registry_loads_expected_reader_class():
    reader = reader_registry["SCF"]

    assert reader.__name__ == "ScfOutputReader"
    assert reader.__module__ == "ase2sprkkr.outputs.readers.scf"


def test_reader_registry_is_lazy_for_discovered_modules(tmp_path, monkeypatch):
    package = _reader_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_reader_registry_lazy",
        {
            "alpha": (
                "from ase2sprkkr.outputs.task_result import KkrOutputReader\n"
                "class AlphaOutputReader(KkrOutputReader):\n"
                "    pass\n"
            )
        },
    )
    registry = OutputReaderRegistry(package)
    module_name = f"{package.__name__}.alpha"

    assert "ALPHA" in registry
    assert module_name not in sys.modules

    reader = registry["ALPHA"]

    assert module_name in sys.modules
    assert reader.__name__ == "AlphaOutputReader"


def test_malformed_reader_module_does_not_fall_back(tmp_path, monkeypatch):
    package = _reader_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_reader_registry_malformed",
        {"broken": "class SomethingElse:\n    pass\n"},
    )
    registry = OutputReaderRegistry(package)

    with pytest.raises(RegistryDefinitionError, match="BrokenOutputReader"):
        registry.class_for_task("BROKEN")


def test_reader_registry_rejects_wrong_export_type(tmp_path, monkeypatch):
    package = _reader_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_reader_registry_wrong_type",
        {"broken": "BrokenOutputReader = object()\n"},
    )
    registry = OutputReaderRegistry(package)

    with pytest.raises(RegistryDefinitionError, match="KkrOutputReader subclass"):
        registry.class_for_task("BROKEN")


def test_reader_dependency_errors_remain_visible(tmp_path, monkeypatch):
    missing_dependency = "_ase2sprkkr_missing_reader_dependency"
    package = _reader_package(
        tmp_path,
        monkeypatch,
        "_ase2sprkkr_reader_registry_dependency",
        {"broken": f"import {missing_dependency}\n"},
    )
    registry = OutputReaderRegistry(package)

    with pytest.raises(ModuleNotFoundError) as exc_info:
        registry.class_for_task("BROKEN")

    assert exc_info.value.name == missing_dependency


def test_reader_selection_facades_delegate_to_registry(monkeypatch):
    expected = object()

    class RecordingRegistry:
        def __init__(self):
            self.tasks = []

        def class_for_task(self, task):
            self.tasks.append(task)
            return expected

    registry = RecordingRegistry()
    monkeypatch.setattr(reader_registry_module, "reader_registry", registry)

    task_definition = TaskDefinition(SimpleNamespace(name="FAKE"))

    assert KkrOutputReader.class_for_task("OUTPUT") is expected
    assert task_definition.output_reader_class() is expected
    assert registry.tasks == ["OUTPUT", "FAKE"]
