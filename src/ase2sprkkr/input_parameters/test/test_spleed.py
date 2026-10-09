"""Dedicated SPLEED energy semantics and reader constraints."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("SPLEED")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "SPLEED"
    assert q.ENERGY.EMINEV() == 11 and q.ENERGY.EMAXEV() == 60
    assert q.ENERGY.NE().tolist() == [200]
    assert q.SPEC_EL.THETA().tolist() == [45., 45.]
    assert q.SPEC_EL.PHI().tolist() == [270., 270.]
    assert q.SPEC_EL.SPOL() == 1
    assert "SPEC_PH" not in text and "EMIN=" not in text
    assert q.SPEC_STR.N_LAYER() == 1
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrspec" and p._definition.mpi is False


@pytest.mark.parametrize("option,value", [
    ("ENERGY.EMINEV", -1.), ("ENERGY.EMAXEV", 5.),
    ("ENERGY.EMIN", 1.), ("ENERGY.NE", [1]),
    ("ENERGY.EWORK_EV", -1.), ("SPEC_EL.NT", 2),
    ("SPEC_EL.THETA", [40., 50.]), ("SPEC_EL.K1", [1., 0.]),
    ("SPEC_EL.POL0", [0, 0, 1]),
])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")


def test_crystal_vector_flag_roundtrip():
    p = parameters()
    p.TASK.CRYS_VECS.set(True)
    text = p.to_string()
    assert "CRYS_VECS" in text
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.CRYS_VECS() is True
