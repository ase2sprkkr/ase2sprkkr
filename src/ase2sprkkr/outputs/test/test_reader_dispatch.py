import importlib
from types import SimpleNamespace

import pytest

from .. import readers
from ..reader_registry import OutputReaderRegistry, reader_registry
from ..task_result import DefaultOutputReader, KkrOutputReader


def _reset_registry_for(monkeypatch, task):
    monkeypatch.setattr(reader_registry, "_module_names", {task.upper(): task})
    monkeypatch.setattr(reader_registry, "_modules", {})
    monkeypatch.setattr(reader_registry, "_values", {})


def test_missing_task_reader_uses_default_reader():
    assert KkrOutputReader.class_for_task("definitely_missing_reader") is DefaultOutputReader


def test_missing_reader_dependency_is_not_hidden(monkeypatch):
    task = "broken"
    missing_dependency = "missing_reader_dependency"
    _reset_registry_for(monkeypatch, task)

    def import_module(name, package=None):
        assert name == f".{task}"
        assert package == readers.__name__
        raise ModuleNotFoundError(
            f"No module named {missing_dependency!r}",
            name=missing_dependency,
        )

    monkeypatch.setattr(importlib, "import_module", import_module)

    with pytest.raises(ModuleNotFoundError) as exc_info:
        KkrOutputReader.class_for_task(task)

    assert exc_info.value.name == missing_dependency


def test_reader_programming_key_error_is_not_treated_as_missing(monkeypatch):
    task = "broken"
    _reset_registry_for(monkeypatch, task)
    monkeypatch.setattr(
        importlib,
        "import_module",
        lambda name, package=None: SimpleNamespace(),
    )

    def broken_value_from_module(self, name, module):
        raise KeyError("reader bug")

    monkeypatch.setattr(
        OutputReaderRegistry,
        "value_from_module",
        broken_value_from_module,
    )

    with pytest.raises(KeyError, match="reader bug"):
        KkrOutputReader.class_for_task(task)
