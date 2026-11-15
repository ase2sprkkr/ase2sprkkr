"""Tests for lightweight and optional atomic-data providers."""

import builtins

import pytest

from ase2sprkkr.sprkkr.atomic_types import AtomicType


def _block_mendeleev_import(monkeypatch):
    real_import = builtins.__import__

    def without_mendeleev(name, *args, **kwargs):
        if name == "mendeleev":
            raise ModuleNotFoundError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_mendeleev)


def test_symbol_number_and_valence_do_not_require_mendeleev(monkeypatch):
    _block_mendeleev_import(monkeypatch)

    nickel = AtomicType("Ni")
    assert nickel.atomic_number == 28
    assert nickel.n_valence == 10
    assert AtomicType(28).symbol == "Ni"


def test_mendeleev_property_remains_optional(monkeypatch):
    _block_mendeleev_import(monkeypatch)
    from ase2sprkkr.physics.atomic_data import mendeleev_element

    mendeleev_element.cache_clear()
    with pytest.raises(ImportError, match=r"ase2sprkkr\[mendeleev\]"):
        _ = AtomicType("Ni").mendeleev


def test_mendeleev_provider_rejects_vacuum():
    from ase2sprkkr.physics.atomic_data import mendeleev_element

    with pytest.raises(ValueError, match="no element"):
        mendeleev_element(0)


def test_mendeleev_property_shares_provider(monkeypatch):
    from ase2sprkkr.physics import atomic_data

    sentinel = object()
    monkeypatch.setattr(atomic_data, "mendeleev_element", lambda z: sentinel if z == 28 else None)
    assert AtomicType("Ni").mendeleev is sentinel
