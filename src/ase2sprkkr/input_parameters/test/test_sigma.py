"""SIGMA grammar and source-derived constraints; no solver execution."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("SIGMA")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip_pickle():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.TASK() == "SIGMA"
    assert q.ENERGY.GRID().tolist() == ["0"] and q.ENERGY.NE().tolist() == [1]
    assert "EMAX=" not in text and "SIGSPROJ" not in text
    assert "SP_NO_ALPHA" not in text and "FULLPOT" not in text
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrchi" and p._definition.mpi is False


def test_sea_projectors_frame_roundtrip():
    p = parameters()
    p.ENERGY.GRID.set([5])
    p.ENERGY.NE.set([30])
    p.TASK.SIGSPROJ.set([1, 2, 3])
    p.TASK.FRAMEGAM.set(45.)
    p.TASK.SP_NO_ALPHA.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.ENERGY.NE().tolist() == [30]
    assert q.TASK.SIGSPROJ().tolist() == [1, 2, 3]
    assert q.TASK.SP_NO_ALPHA() is True and q.TASK.FRAMEGAM() == 45.


def test_nabla_and_nonmag():
    p = parameters()
    p.CONTROL.NONMAG.set(True)
    p.TASK.SIGSPROJ.set([9])
    p.validate("save")


@pytest.mark.parametrize("option,value,match", [
    ("TASK.SIGSPROJ", [2, 1], "first"),
    ("TASK.SIGSPROJ", [1, 1], "duplicates"),
    ("TASK.SIGSPROJ", [1, 17], None),
    ("ENERGY.NE", [30], "GRID=5"),
    ("ENERGY.GRID", [5], "GRID=5"),
    ("ENERGY.SPLITSS", True, "SPLITSS"),
    ("ENERGY.ImE", 0., "ImE"),
    ("ENERGY.ImE", float("nan"), "ImE"),
    ("SCF.FULLPOT", True, "ASA"),
    ("MODE.MODE", "SREL", "FREL"),
    ("TAU.BZINT", "TETRA", "POINTS"),
])
def test_invalid(option, value, match):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError), match=match):
        p.set({option: value})
        p.validate("save")
