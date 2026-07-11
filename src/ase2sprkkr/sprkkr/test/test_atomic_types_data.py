"""Tests for lightweight and optional atomic-data providers."""

import builtins
import pytest

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


from ase2sprkkr.sprkkr.atomic_types import AtomicType


def test_symbol_number_and_valence_do_not_require_mendeleev(monkeypatch):
    real_import = builtins.__import__

    def without_mendeleev(name, *args, **kwargs):
        if name == "mendeleev":
            raise ModuleNotFoundError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_mendeleev)

    nickel = AtomicType("Ni")
    assert nickel.atomic_number == 28
    assert nickel.n_valence == 10
    assert AtomicType(28).symbol == "Ni"


def test_mendeleev_property_remains_optional(monkeypatch):
    real_import = builtins.__import__

    def without_mendeleev(name, *args, **kwargs):
        if name == "mendeleev":
            raise ModuleNotFoundError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_mendeleev)
    AtomicType._mendeleev_module = None

    with pytest.raises(ImportError, match=r"ase2sprkkr\[mendeleev\]"):
        _ = AtomicType("Ni").mendeleev
