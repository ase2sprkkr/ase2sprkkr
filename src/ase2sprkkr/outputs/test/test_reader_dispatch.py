from types import SimpleNamespace

import pytest

from ...common.module_registry import (
    RegistryDefinitionError,
    RegistryEntryNotFound,
)
from .. import reader_registry as reader_registry_module
from ..reader_registry import reader_registry
from ..readers.scf import ScfOutputReader
from ..task_result import DefaultOutputReader, KkrOutputReader


class _MissingRegistry:
    def __getitem__(self, name):
        raise RegistryEntryNotFound(name)


class _BrokenRegistry:
    def __getitem__(self, name):
        raise KeyError("reader bug")


def test_missing_task_reader_uses_default_reader():
    assert KkrOutputReader.class_for_task("definitely_missing_reader") is DefaultOutputReader


def test_reader_registry_loads_expected_reader_class():
    assert reader_registry["SCF"] is ScfOutputReader


def test_reader_module_without_expected_class_is_malformed():
    module = SimpleNamespace(__name__="fake_reader")

    with pytest.raises(RegistryDefinitionError, match="BrokenOutputReader"):
        reader_registry.value_from_module("BROKEN", module)


def test_reader_registry_rejects_wrong_export_type():
    module = SimpleNamespace(__name__="fake_reader")

    with pytest.raises(RegistryDefinitionError, match="KkrOutputReader subclass"):
        reader_registry.prepare_value("BROKEN", module, object())


def test_class_for_task_falls_back_only_for_missing_entry(monkeypatch):
    monkeypatch.setattr(
        reader_registry_module,
        "reader_registry",
        _MissingRegistry(),
    )

    assert KkrOutputReader.class_for_task("missing") is DefaultOutputReader


def test_class_for_task_keeps_reader_programming_errors_visible(monkeypatch):
    monkeypatch.setattr(
        reader_registry_module,
        "reader_registry",
        _BrokenRegistry(),
    )

    with pytest.raises(KeyError, match="reader bug"):
        KkrOutputReader.class_for_task("broken")
