"""Tests for selection of the empty-spheres backend."""

import warnings

import pytest

from .. import empty_spheres as module
from ..xband import spheres


def test_unavailable_es_finder_warns_only_for_explicit_selection(monkeypatch):
    expected = object()
    monkeypatch.setattr(module.es_finder, "is_enabled", False)
    monkeypatch.setattr(spheres, "empty_spheres", lambda atoms, **kwargs: expected)

    with pytest.warns(RuntimeWarning, match="Falling back.*inhouse"):
        result = module.empty_spheres(object(), method="es_finder")

    assert result is expected

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        automatic = module.empty_spheres(object(), method="auto")
    assert automatic is expected
    assert caught == []


def test_available_es_finder_is_selected(monkeypatch):
    expected = object()
    monkeypatch.setattr(module.es_finder, "is_enabled", True)
    monkeypatch.setattr(module.es_finder, "empty_spheres", lambda atoms, **kwargs: expected)

    assert module.empty_spheres(object(), method="es_finder") is expected


