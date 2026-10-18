"""Input checks only; no numerical positron preparation."""
import io
import pickle
import warnings
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors, DataValidityWarning


def parameters():
    p = InputParameters.create("POSANITAU")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.ENERGY.EMIN.set(-0.1)
    return p


def test_defaults_and_discovery(tmp_path):
    p = parameters()
    with pytest.warns(DataValidityWarning, match="band") as caught:
        text = p.to_string()
    assert "POSVC" in str(caught[0].message)
    assert "GRID={5}" in text and "NE={30}" in text
    assert "POSVC=BM" in text and "EMIN=-0.1" in text
    assert "EMAX=" not in text and "SPLITSS" not in text
    filename = tmp_path / "prep.in"
    filename.write_text(text)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DataValidityWarning)
        q = InputParameters.from_file(filename)
        assert q.TASK.TASK() == "POSANITAU"
        assert pickle.loads(pickle.dumps(q)).to_string() == text
    assert q._definition.mpi is False


@pytest.mark.parametrize("posvc", ["BM", "DR", "IPM"])
def test_posvc(posvc):
    p = parameters()
    p.TASK.POSVC.set(posvc)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DataValidityWarning)
        q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.TASK.POSVC() == posvc


@pytest.mark.parametrize("section,name,value", [
    ("MODE", "MODE", "SREL"), ("MODE", "MODE", "SP-SREL"), ("MODE", "MODE", "NREL"),
    ("SCF", "FULLPOT", True), ("ENERGY", "GRID", [0]), ("ENERGY", "NE", [30,30]),
    ("ENERGY", "EMINEV", 1.), ("ENERGY", "EMAXEV", 1.), ("TAU", "BZINT", "CLUSTER"),
])
def test_constraints(section,name,value):
    p = parameters()
    p[section][name].set(value)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DataValidityWarning)
        with pytest.raises((DataValidityError, DataValidityErrors)):
            p.validate("save")


def test_nonmag_allowed():
    p = parameters()
    p.CONTROL.NONMAG.set(True)
    with pytest.warns(DataValidityWarning):
        p.validate("save")


def test_upper_bound():
    p = parameters()
    p.ENERGY.EMAX.set(0.2)
    with pytest.warns(DataValidityWarning):
        text = p.to_string()
    assert "EMAX=0.2" in text
