"""DOS input regressions; no numerical solver execution."""
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors, DataValidityWarning

ERRORS = (DataValidityError, DataValidityErrors)


def parameters():
    p = InputParameters.create("DOS")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_roundtrip(tmp_path):
    p = parameters()
    text = p.to_string()
    assert "GRID={3}" in text and "NE={300}" in text
    assert "EMAX=1.0" in text and "ImE=0.01" in text
    filename = tmp_path / "dos.in"
    filename.write_text(text)
    q = InputParameters.from_file(filename)
    assert q.TASK.TASK() == "DOS"
    assert pickle.loads(pickle.dumps(q)).to_string() == text
    assert q._definition.mpi is True


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_modes_roundtrip(mode):
    p = parameters()
    p.MODE.MODE.set(mode)
    q = p._definition.read_from_string(p.to_string())
    assert q.MODE.MODE() == mode


@pytest.mark.parametrize("grid", [2,3])
def test_single_point_unsafe(grid):
    p = parameters()
    p.ENERGY.GRID.set([grid])
    p.ENERGY.NE.set([1])
    with pytest.raises(ERRORS, match="two energy"):
        p.validate("save")


def test_grid_zero_single_point():
    p = parameters()
    p.ENERGY.GRID.set([0])
    p.ENERGY.NE.set([1])
    p.validate("save")


@pytest.mark.parametrize("ne", [0,-1])
def test_nonpositive_ne(ne):
    p = parameters()
    with pytest.raises(ERRORS):
        p.ENERGY.NE.set([ne])


def test_zero_ime_warning():
    p = parameters()
    p.ENERGY.ImE.set(0.)
    with pytest.warns(DataValidityWarning, match="real axis"):
        p.validate("save")


def test_negative_ime():
    p = parameters()
    p.ENERGY.ImE.set(-0.01)
    with pytest.raises(ERRORS, match="nonnegative"):
        p.validate("save")


def test_relative_bounds():
    p = parameters()
    p.ENERGY.EMINEV.set(-10.)
    p.ENERGY.EMAXEV.set(5.)
    text = p.to_string()
    assert "EMINEV=-10.0" in text and "EMAXEV=5.0" in text
    assert "EMIN=" not in text and "EMAX=" not in text
    p._definition.read_from_string(text).validate("save")
