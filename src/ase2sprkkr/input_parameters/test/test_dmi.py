"""DMI preset and native DMIVECTOR constraints; no numerical solver run."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("DMI")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_defaults_and_roundtrip(tmp_path):
    p = parameters()
    assert p.name == "jxc"
    assert p.TASK.TASK() == "JXC"
    assert not p.TASK.JXC() and p.TASK.DMI()
    text = p.to_string()
    assert "\n\tDMI\n" in text and "\n\tJXC\n" not in text
    assert "CLURAD=2.2" in text
    assert "NE={30}" in text and "ImE=0.01" in text
    assert "EMAX=" not in text and "QREPONLY" not in text
    filename = tmp_path / "dmi.in"
    filename.write_text(text)
    q = InputParameters.from_file(filename)
    assert q.TASK.TASK() == "JXC"
    assert not q.TASK.JXC() and q.TASK.DMI()
    assert q.ENERGY.GRID().tolist() == ["5"]
    assert q.ENERGY.NE().tolist() == [30]
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrgen"
    assert p._definition.mpi is True


def test_copy_preserves_dmi_preset():
    p = parameters().copy(copy_values=True)
    assert not p.TASK.JXC() and p.TASK.DMI()
    assert p.ENERGY.NE().tolist() == [30]
    assert p.ENERGY.ImE() == 0.01


def test_select_shells_and_qreponly_roundtrip():
    p = parameters()
    p.TASK.NSHLCLU.set(0)
    p.TASK.SELECT.set(["Fe", "Co"])
    p.TASK.QREPONLY.set(True)
    text = p.to_string()
    assert "NSHLCLU=0" in text and "QREPONLY" in text
    q = p._definition.read_from_file(io.StringIO(text))
    assert q.TASK.SELECT().tolist() == ["Fe", "Co"]
    assert q.TASK.NSHLCLU() == 0


def test_exclude():
    p = parameters()
    p.TASK.EXCLUDE.set(["Cu"])
    p.validate("save")
    p.TASK.SELECT.set(["Fe"])
    with pytest.raises((DataValidityError, DataValidityErrors), match="mutually exclusive"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_reject_non_frel(mode):
    p = parameters()
    p.MODE.MODE.set(mode)
    with pytest.raises((DataValidityError, DataValidityErrors), match="magnetic FREL"):
        p.validate("save")


def test_reject_nonmag():
    p = parameters()
    p.CONTROL.NONMAG.set(True)
    with pytest.raises((DataValidityError, DataValidityErrors), match="magnetic"):
        p.validate("save")


@pytest.mark.parametrize("grid", [[3], [5, 5]])
def test_reject_contours(grid):
    p = parameters()
    p.ENERGY.GRID.set(grid)
    with pytest.raises((DataValidityError, DataValidityErrors)):
        p.validate("save")


def test_reject_splitss():
    p = parameters()
    p.ENERGY.GRID.set([5, 5])
    p.ENERGY.NE.set([30, 30])
    p.ENERGY.SPLITSS.set(True)
    with pytest.raises((DataValidityError, DataValidityErrors), match="one GRID=5"):
        p.validate("save")


@pytest.mark.parametrize("bzint", ["WEYL", "CLUSTER"])
def test_reject_other_workflows(bzint):
    p = parameters()
    p.TAU.BZINT.set(bzint)
    with pytest.raises((DataValidityError, DataValidityErrors), match="BZINT=POINTS"):
        p.validate("save")


@pytest.mark.parametrize("ne", [0, -1])
def test_reject_invalid_ne(ne):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors)):
        p.ENERGY.NE.set([ne])
        p.validate("save")
