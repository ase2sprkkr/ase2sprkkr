import importlib

import pytest

from .. import readers
from ..task_result import DefaultOutputReader, KkrOutputReader


def test_missing_task_reader_uses_default_reader():
    assert KkrOutputReader.class_for_task("definitely_missing_reader") is DefaultOutputReader


def test_missing_reader_dependency_is_not_hidden(monkeypatch):
    task = "broken"
    missing_dependency = "missing_reader_dependency"

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
