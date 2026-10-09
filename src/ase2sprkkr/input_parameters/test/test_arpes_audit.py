"""ARPES native reader domains, retaining valid single-point scans."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("ARPES")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.ENERGY.EWORK_EV() == 4.2
    assert q.SPEC_STR.N_LAYDBL().tolist() == [10, 10]
    assert q.SPEC_EL.NT() == 1 and q.SPEC_EL.NP() == 1
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrspec"


def test_single_k_point_roundtrip():
    p = parameters()
    p.SPEC_EL.K1.set([0., 0.])
    p.SPEC_EL.NK1.set(1)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.SPEC_EL.NK1() == 1


def test_missing_k_count():
    p = parameters()
    p.SPEC_EL.K1.set([1., 0.])
    with pytest.raises((DataValidityError, DataValidityErrors), match="requires NK1"):
        p.validate("save")


@pytest.mark.parametrize("option,value", [
    ("SPEC_EL.NT", 0), ("SPEC_EL.NP", -1),
    ("SPEC_EL.NK1", 0), ("SPEC_STR.N_LAYER", 0),
    ("SPEC_STR.NLAT_G_VEC", -1), ("SPEC_STR.N_LAYDBL", [1]),
    ("ENERGY.EWORK_EV", -1.), ("ENERGY.IMV_INI_EV", float("nan")),
    ("SPEC_PH.EPHOT", 0.), ("TASK.INPVER", 2),
])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")


def test_ignored_extra_k_directions():
    p = parameters()
    p.SPEC_EL.K3.set([1., 0.])
    with pytest.raises((DataValidityError, DataValidityErrors), match="does not read"):
        p.validate("save")
