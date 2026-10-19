"""Input regression checks; no numerical COMPTON solver execution."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors, DataValidityWarning

ERRORS = (DataValidityError, DataValidityErrors)


def parameters():
    p = InputParameters.create("COMPTON")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.TASK.PNMAX.set(3.0)
    return p


def test_defaults_and_discovery(tmp_path):
    p = parameters()
    text = p.to_string()
    assert text.count("\tSPLITSS\n") == 2
    assert "GRID={5,5}" in text and "NE={30,50}" in text
    assert "EMAX=" not in text and "ImE=0.01" in text
    filename = tmp_path / "compton.in"
    filename.write_text(text)
    q = InputParameters.from_file(filename)
    assert q.TASK.TASK() == "COMPTON"
    assert q.TASK.NPMAX() == 250
    assert pickle.loads(pickle.dumps(q)).to_string() == text
    assert q._definition.mpi is True


def test_default_zero_warning():
    p = parameters()
    p.TASK.PNMAX.clear()
    with pytest.warns(DataValidityWarning, match="one point"):
        p.validate("save")


def test_vectors_and_paths_roundtrip():
    p = parameters()
    p.TASK.PNVEC.set([0., 0., 1.])
    p.TASK.PPVEC1.set([1., 0., 0.])
    p.TAU.TAUDIR.set("Fe.tau.dir")
    p.TAU.TAUPATH.set("cache")
    p.TASK.MEMODE.set(2)
    p.CONTROL.NONMAG.set(True)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.TASK.PNVEC().tolist() == [0., 0., 1.]
    assert q.TAU.TAUPATH() == "cache"
    assert q.CONTROL.NONMAG()


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_mode(mode):
    p = parameters()
    p.MODE.MODE.set(mode)
    with pytest.raises(ERRORS, match="FREL"):
        p.validate("save")


@pytest.mark.parametrize("section", ["ENERGY", "TASK"])
def test_splitss(section):
    p = parameters()
    p[section].SPLITSS.set(False)
    with pytest.raises(ERRORS, match="SPLITSS"):
        p.validate("save")


@pytest.mark.parametrize("values", [dict(PNVEC=[0.,0.,0.]), dict(PPVEC1=[1.,0.,0.]),
                                     dict(PNVEC=[0.,0.,1.], PPVEC1=[1.,0.,1.]),
                                     dict(PNVEC=[0.,0.,1.], PPVEC2=[0.,1.,0.])])
def test_invalid_vectors(values):
    p = parameters()
    for name, value in values.items():
        p.TASK[name].set(value)
    with pytest.raises(ERRORS):
        p.validate("save")


def test_transverse_pair():
    p = parameters()
    p.TASK.PPVEC1.set([1.,0.,0.])
    p.TASK.PPVEC2.set([0.,1.,0.])
    p.validate("save")


def test_shifted_grid():
    p = parameters()
    p.TASK.SHFT.set(True)
    with pytest.raises(ERRORS, match="odd"):
        p.validate("save")
    p.TASK.NPN.set(51)
    p.TASK.NPP.set(51)
    p.validate("save")
    p.TASK.SHFT_PP.set(True)
    with pytest.raises(ERRORS, match="at most one"):
        p.validate("save")


def test_emd2d():
    p = parameters()
    p.TASK.EMD2D.set(True)
    with pytest.raises(ERRORS, match="scaling twice"):
        p.validate("save")
    p.TASK.PNMAX.set(0.)
    p.validate("save")


@pytest.mark.parametrize("name,value", [("NPMAX",1), ("NPMAX",10001), ("MEMODE",3), ("NPN",0), ("USE_PEQUIV",3)])
def test_invalid_counts(name,value):
    p = parameters()
    with pytest.raises(ERRORS):
        p.TASK[name].set(value)


def test_energy_ne():
    p = parameters()
    with pytest.raises(ERRORS):
        p.ENERGY.NE.set([0,50])


def test_contour_and_bz():
    p = parameters()
    p.ENERGY.GRID.set([3,5])
    with pytest.raises(ERRORS, match="GRID=5"):
        p.validate("save")
    p.ENERGY.GRID.set([5,5])
    p.TAU.BZINT.set("CLUSTER")
    with pytest.raises(ERRORS, match="POINTS"):
        p.validate("save")
