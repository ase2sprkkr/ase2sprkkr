"""GILBERT input constraints, without numerical response calculation."""
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors

ERRORS = (DataValidityError, DataValidityErrors)


def parameters():
    p = InputParameters.create("GILBERT")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_and_discovery(tmp_path):
    p = parameters()
    text = p.to_string()
    assert "GRID={0}" in text and "NE={1}" in text
    assert "EMIN=" not in text and "EMAX=" not in text
    assert "SETFLUCT" not in text and "GILBERT" in text
    path = tmp_path / "gil.in"
    path.write_text(text)
    q = InputParameters.from_file(path)
    assert q.TASK.TASK() == "GILBERT"
    assert pickle.loads(pickle.dumps(q)).to_string() == text
    assert q.task_definition.executable == "kkrchi"


def test_mlin_spin():
    p = parameters()
    p.TASK.SETFLUCT.set("MLIN")
    p.TASK.TMPMIN.set(0.)
    p.TASK.TMPMAX.set(300.)
    p.TASK.NTMP.set(4)
    p.TASK.NFPHI.set(4)
    p.TASK.NFTET.set(4)
    with pytest.raises(ERRORS, match="DMAGDT"):
        p.validate("save")
    p.TASK.DMAGDT.set(0.001)
    p.TASK.TDEBYE.set(400.)
    p.TASK.NVIBRA.set(14)
    p.TASK.INITSVEC.set(0)
    q = p._definition.read_from_string(p.to_string())
    assert q.TASK.DMAGDT() == 0.001
    assert q.TASK.TDEBYE() == 400.


def test_m_t_file():
    p = parameters()
    p.TASK.SETFLUCT.set("M_T")
    with pytest.raises(ERRORS, match="FLUCTFIL"):
        p.validate("save")
    p.TASK.FLUCTFIL.set("mag.dat")
    p.validate("save")
    p.TASK.NTMP.set(4)
    with pytest.raises(ERRORS, match="temperature mesh"):
        p.validate("save")


def test_lattice_complete_mesh():
    p = parameters()
    p.TASK.TMPMIN.set(0.)
    with pytest.raises(ERRORS, match="together"):
        p.validate("save")
    p.TASK.TMPMAX.set(300.)
    p.TASK.NTMP.set(4)
    p.validate("save")


@pytest.mark.parametrize("name", ["NTMP","NVIBRA","NFTET","NFPHI"])
def test_positive_counts(name):
    p = parameters()
    with pytest.raises(ERRORS):
        p.TASK[name].set(0)


@pytest.mark.parametrize("mode", ["NREL","SREL","SP-SREL"])
def test_magnetic_frel(mode):
    p = parameters()
    p.MODE.MODE.set(mode)
    with pytest.raises(ERRORS, match="FREL"):
        p.validate("save")


def test_nonmag():
    p = parameters()
    p.CONTROL.NONMAG.set(True)
    with pytest.raises(ERRORS, match="FREL"):
        p.validate("save")


def test_wrong_grid():
    p = parameters()
    p.ENERGY.GRID.set([3])
    with pytest.raises(ERRORS, match="GRID=0"):
        p.validate("save")


def test_missing_mlin_temperature():
    p = parameters()
    p.TASK.SETFLUCT.set("MLIN")
    with pytest.raises(ERRORS, match="TMPMIN"):
        p.validate("save")
