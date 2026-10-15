"""ARPES absolute/relative bounds follow INIT_MOD_ENERGY and SPEC_INPUT."""
from io import StringIO
import pickle

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.common.warnings import DataValidityError


def parameters():
    result = InputParameters.create("arpes")
    result.CONTROL.POTFIL = "Fe.pot"
    return result


def test_relative_defaults():
    ip = parameters()
    assert ip.ENERGY.EMIN() is None
    assert ip.ENERGY.EMAX() is None
    assert ip.ENERGY.EMINEV() == -8.
    assert ip.ENERGY.EMAXEV() == 5.
    ip.validate("save")


def test_absolute_bounds_suppress_relative_defaults_and_roundtrip():
    ip = parameters()
    ip.ENERGY.set({"EMIN": .2, "EMAX": .8})
    assert ip.ENERGY.EMINEV() is None
    assert ip.ENERGY.EMAXEV() is None
    source = ip.to_string(validate=True)
    assert "EMINEV=" not in source
    assert "EMAXEV=" not in source
    restored = parameters()
    restored.read_from_file(StringIO(source))
    for candidate in (restored, ip.copy(), pickle.loads(pickle.dumps(ip))):
        assert candidate.ENERGY.EMIN() == .2
        assert candidate.ENERGY.EMAX() == .8
        assert candidate.ENERGY.EMINEV() is None
        assert candidate.ENERGY.EMAXEV() is None
        candidate.validate("save")
    ip.ENERGY.set({"EMIN": None, "EMAX": None})
    assert ip.ENERGY.EMINEV() == -8.
    assert ip.ENERGY.EMAXEV() == 5.


@pytest.mark.parametrize("name", ["EMIN", "EMAX"])
def test_mixed_reference_is_editable_but_not_saveable(name):
    ip = parameters()
    ip.ENERGY[name].set(.2)
    with pytest.raises(DataValidityError, match="must be supplied together"):
        ip.validate("save")


def test_explicit_relative_pair_overrides_absolute_as_in_sprkkr():
    ip = parameters()
    ip.ENERGY.set({"EMIN": .2, "EMAX": .8, "EMINEV": -3., "EMAXEV": 2.})
    assert ip.ENERGY.EMINEV() == -3.
    assert ip.ENERGY.EMAXEV() == 2.
    ip.validate("save")
