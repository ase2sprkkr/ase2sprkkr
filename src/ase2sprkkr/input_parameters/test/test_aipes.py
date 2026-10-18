"""Dedicated AIPES definition and native workflow scope."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("AIPES")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "AIPES"
    assert q.ENERGY.EMINEV() == -8 and q.ENERGY.EMAXEV() == 0
    assert q.SPEC_PH.EPHOT() == 25
    assert q.SPEC_EL.THETA().tolist() == [0., 0.]
    assert q.SPEC_EL.SPOL() == 2
    assert "CRYS_VEC\n" not in text
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p.task_definition.executable == "kkrspec" and p.task_definition.mpi is False


def test_work_function_and_photons():
    p = parameters()
    p.ENERGY.EWORK_EV.set(5.)
    p.SPEC_PH.EPHOT.set(1253.6)
    p.TASK.CRYS_VECS.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.ENERGY.EWORK_EV() == 5 and q.SPEC_PH.EPHOT() == 1253.6
    assert q.TASK.CRYS_VECS() is True


@pytest.mark.parametrize("option,value", [
    ("SPEC_EL.NT", 2), ("SPEC_EL.NP", 0),
    ("SPEC_EL.THETA", [0., 45.]), ("SPEC_EL.K1", [1., 0.]),
    ("ENERGY.EWORK_EV", -1.), ("ENERGY.IMV_FIN_EV", float("nan")),
    ("SPEC_PH.EPHOT", 0.), ("ENERGY.SPLITSS", True), ("ENERGY.NE", [0]),
])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")
