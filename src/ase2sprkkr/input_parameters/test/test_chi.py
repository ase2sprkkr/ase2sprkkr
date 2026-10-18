"""CHI grammar and native restrictions; no numerical SPR-KKR execution."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("CHI")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "CHI"
    assert q.ENERGY.NE().tolist() == [30]
    assert q.ENERGY.GRID().tolist() == ["5"]
    assert q.MODE.NCHIITER() == 1 and "NCHIITER=1" in text
    assert "EMAX=" not in text and "NOSYMCHI" not in text
    assert "FULLPOT" not in text
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrchi" and p._definition.mpi is False


def test_mode_controls_roundtrip():
    p = parameters()
    p.MODE.PRINT.set(2)
    p.TASK.CHIPRINT.set(3)
    p.TASK.NOSYMCHI.set(True)
    p.MODE.OP.set("BROOKS")
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.MODE.PRINT() == 2
    assert q.TASK.CHIPRINT() == 3
    assert q.TASK.NOSYMCHI() is True
    assert q.MODE.OP() == "BROOKS"


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
    ("MODE.OP", "DMFT", "OP=NONE"),
    ("MODE.NCHIITER", 2, None),
    ("ENERGY.NE", [0], None),
])
def test_invalid(option, value, match):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError), match=match):
        p.set({option: value})
        p.validate("save")
