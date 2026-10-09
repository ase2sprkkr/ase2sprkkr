"""Compatibility warning for legacy KKRSPEC aliases."""
import warnings
import pytest
from ..input_parameters import InputParameters

@pytest.mark.parametrize("task", ["AIPES", "BAND", "SPLEED"])
def test_legacy_alias_warns(task):
    p = InputParameters.create("ARPES")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.TASK.TASK.set(task)
    with pytest.warns(DeprecationWarning, match="task-specific defaults"):
        p.validate("save")

def test_arpes_does_not_warn():
    p = InputParameters.create("ARPES")
    p.CONTROL.POTFIL.set("Fe.pot")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p.validate("save")
    assert not [w for w in caught if issubclass(w.category, DeprecationWarning)]
