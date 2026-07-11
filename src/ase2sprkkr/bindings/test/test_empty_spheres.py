"""Tests for selection of the empty-spheres backend."""

import pytest

from .. import empty_spheres as module


def test_requested_unavailable_es_finder_falls_back(monkeypatch):
    expected = object()
    monkeypatch.setattr(module.es_finder, "is_enabled", False)
    monkeypatch.setattr(module.spheres, "empty_spheres", lambda atoms, **kwargs: expected)

    with pytest.warns(RuntimeWarning, match="Falling back.*inhouse"):
        result = module.empty_spheres(object(), method="es_finder")

    assert result is expected


def test_available_es_finder_is_selected(monkeypatch):
    expected = object()
    monkeypatch.setattr(module.es_finder, "is_enabled", True)
    monkeypatch.setattr(module.es_finder, "empty_spheres", lambda atoms, **kwargs: expected)

    assert module.empty_spheres(object(), method="es_finder") is expected


def test_auto_silently_uses_inhouse_when_es_finder_is_unavailable(monkeypatch):
    expected = object()
    monkeypatch.setattr(module.es_finder, "is_enabled", False)
    monkeypatch.setattr(module.spheres, "empty_spheres", lambda atoms, **kwargs: expected)

    assert module.empty_spheres(object(), method="auto") is expected
