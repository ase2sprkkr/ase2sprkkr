"""KKRSPEC BAND grammar and native root-filter safeguards."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("BAND")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "BAND"
    assert q.ENERGY.IMV_INI_EV() == 0 and q.ENERGY.IMV_FIN_EV() == 0
    assert q.SPEC_EL.NK1() == 1
    assert q.SPEC_BAND.EIG_IME() == 0.01 and "EIG_IME=0.01" in text
    assert q.SPEC_BAND.NEIGMAX() == 100000
    assert "CALC_TRANS_PROB" not in text
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p.task_definition.executable == "kkrspec" and p.task_definition.mpi is False


def test_momentum_and_transmission_roundtrip():
    p = parameters()
    p.SPEC_EL.K1.set([1., 0.])
    p.SPEC_EL.NK1.set(20)
    p.SPEC_BAND.CALC_TRANS_PROB.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.SPEC_EL.NK1() == 20 and q.SPEC_BAND.CALC_TRANS_PROB() is True


@pytest.mark.parametrize("option,value", [
    ("SPEC_BAND.EIG_IME", 0.), ("SPEC_BAND.EIG_IME", float("nan")),
    ("SPEC_BAND.NEIGMAX", 0), ("SPEC_EL.NK1", 0),
    ("SPEC_EL.K3", [1., 0.]), ("ENERGY.IMV_INI_EV", -1.),
    ("ENERGY.NE", [0]), ("ENERGY.SPLITSS", True),
])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")
