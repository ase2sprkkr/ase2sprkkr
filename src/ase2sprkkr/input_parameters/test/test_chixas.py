"""CHIXAS external data and grammar; no numerical solver run."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("CHIXAS")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.TASK.NEF.set([[1., 2., 3.]])
    p.TASK.CHISPN.set([1.])
    p.TASK.CHIVV.set([0.2])
    return p


def test_defaults_and_data_roundtrip():
    p = parameters()
    text = p.to_string()
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.ITXRAY() == 1 and q.TASK.CL() == "2P"
    assert "NEF1=" in text and "CHISPN=" in text
    assert q.TASK.NEF(all_values=True)[0].tolist() == [1., 2., 3.]
    assert q.ENERGY.EMAX() == 4 and q.ENERGY.NE().tolist() == [180]
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p.task_definition.executable == "kkrchi" and p.task_definition.mpi is False


def test_missing_external_data():
    p = InputParameters.create("CHIXAS")
    p.CONTROL.POTFIL.set("Fe.pot")
    with pytest.raises((DataValidityError, DataValidityErrors), match="required"):
        p.validate("save")


def test_two_types_roundtrip():
    p = parameters()
    p.TASK.NEF.set([[1., 2.], [3., 4.]])
    p.TASK.CHISPN.set([1., 2.])
    p.TASK.CHIVV.set([0.2, 0.3])
    p.TASK.ITXRAY.set(2)
    p.TASK.ABSUNIT.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.TASK.ITXRAY() == 2 and q.TASK.ABSUNIT() is True
    assert q.TASK.NEF(all_values=True)[1].tolist() == [3., 4.]


@pytest.mark.parametrize("option,value,match", [
    ("TASK.ITXRAY", 2, "type count"),
    ("TASK.CHIVV", [1., 2.], "same contiguous"),
    ("TASK.NEF", [[1.], [2.]], "same contiguous"),
    ("TASK.NEF", [[-1.]], "nonnegative"),
    ("TASK.CHISPN", [float("nan")], "finite"),
    ("ENERGY.NE", [1], None),
    ("ENERGY.GRID", [5], "spectral mesh"),
    ("ENERGY.ImE", 0., "ImE"),
    ("SCF.FULLPOT", True, "FULLPOT"),
    ("MODE.MODE", "SREL", "FREL"),
])
def test_invalid(option, value, match):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError), match=match):
        p.set({option: value})
        p.validate("save")
