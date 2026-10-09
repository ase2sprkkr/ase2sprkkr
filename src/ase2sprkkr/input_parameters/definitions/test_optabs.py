"""OPTABS source-derived constraints; no numerical solver run."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("OPTABS")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "OPTABS"
    assert q.TASK.NEB() == 30 and q.TASK.NOMEGA() == 100
    assert q.TASK.TEMPEL() == 300 and q.TASK.NKT() == 6
    assert "OMEGA1=0.0" in text
    assert "EMIN=" not in text and "EMAX=" not in text
    assert q.ENERGY.GRID().tolist() == ["0"] and q.ENERGY.NE().tolist() == [1]
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrchi" and p._definition.mpi is False


def test_frequency_resolution_roundtrip():
    p = parameters()
    p.TASK.NEB.set(150)
    p.TASK.NOMEGA.set(150)
    p.TASK.FRAMEGAM.set(45.)
    p.CONTROL.NONMAG.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.TASK.NEB() == 150 and q.TASK.FRAMEGAM() == 45.


@pytest.mark.parametrize("option,value,match", [
    ("TASK.TEMPEL", 0., "KBT"),
    ("TASK.TEMPEL", 1., "KBT"),
    ("TASK.OMEGA2", 0., "OMEGA2"),
    ("TASK.OMEGA2", float("nan"), "OMEGA2"),
    ("TASK.OMEGA1", 1., "OMEGA1"),
    ("TASK.NOMEGA", 3, "overflow"),
    ("TASK.NEB", 2, None),
    ("SCF.FULLPOT", True, "FULLPOT"),
    ("MODE.MODE", "SREL", "FREL"),
    ("ENERGY.GRID", [3], "initialization"),
    ("ENERGY.ImE", 0., "ImE"),
])
def test_invalid(option, value, match):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError), match=match):
        p.set({option: value})
        p.validate("save")
