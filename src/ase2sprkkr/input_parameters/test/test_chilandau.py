"""CHILANDAU grammar and source-derived constraints; no solver execution."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("CHILANDAU")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "CHILANDAU"
    assert q.ENERGY.GRID().tolist() == ["5"]
    assert q.ENERGY.NE().tolist() == [30]
    assert "EMAX=" not in text and "CLURAD=2.2" in text
    assert "FULLPOT" not in text and "QREPONLY" not in text
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrchi"
    assert p._definition.mpi is False


def test_shell_selection_roundtrip():
    p = parameters()
    p.TASK.NSHLCLU.set(0)
    p.TASK.SELECT.set(["Fe", "Co"])
    p.TASK.PRINT.set(2)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.TASK.NSHLCLU() == 0
    assert q.TASK.SELECT().tolist() == ["Fe", "Co"]
    assert q.TASK.PRINT() == 2


def test_nonmag_allowed():
    p = parameters()
    p.CONTROL.NONMAG.set(True)
    p.validate("save")


@pytest.mark.parametrize("option,value,match", [
    ("MODE.MODE", "SREL", "FREL"),
    ("SCF.FULLPOT", True, "FULLPOT"),
    ("TAU.BZINT", "TETRA", "POINTS"),
    ("ENERGY.GRID", [3], "GRID=5"),
    ("ENERGY.SPLITSS", True, "GRID=5"),
    ("TASK.CLURAD", -1., "CLURAD"),
    ("TASK.CLURAD", float("nan"), "CLURAD"),
    ("ENERGY.NE", [0], None),
])
def test_invalid(option, value, match):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError), match=match):
        p.set({option: value})
        p.validate("save")


def test_exclusive_selectors():
    p = parameters()
    p.TASK.SELECT.set(["Fe"])
    p.TASK.EXCLUDE.set(["Co"])
    with pytest.raises((DataValidityError, DataValidityErrors), match="mutually exclusive"):
        p.validate("save")
