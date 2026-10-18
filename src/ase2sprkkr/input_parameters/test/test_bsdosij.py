import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("BSDOSIJ")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert p.is_mpi() and p._definition.executable == "kkrgen"
    assert p.ENERGY.GRID().tolist() == ["3"] and p.TAU.BZINT() == "POINTS"
    text = p.to_string(validate=True)
    assert "EMAX=" not in text
    path = tmp_path / "bsdosij.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


@pytest.mark.parametrize("section,name,value", [
    ("ENERGY", "GRID", 5), ("ENERGY", "GRID", [3, 5]),
    ("ENERGY", "NE", [30, 10]), ("ENERGY", "SPLITSS", True),
    ("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"),
    ("MODE", "MODE", "SP-SREL"), ("CONTROL", "NONMAG", True),
    ("TAU", "BZINT", "WEYL"),
])
def test_restrictions(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


def test_selection(tmp_path):
    p = parameters()
    p.TASK.NSHLCLU = 2
    p.TASK.SELECT = ["Fe", "Co"]
    p.TASK.QREPONLY = True
    path = tmp_path / "selected.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.SELECT().tolist() == ["Fe", "Co"]
    assert q.TASK.NSHLCLU() == 2 and q.TASK.QREPONLY()
    p.TASK.EXCLUDE = ["Ni"]
    with pytest.raises(DataValidityError, match="mutually exclusive"):
        p.validate("save")


@pytest.mark.parametrize("value", [0, -1])
def test_invalid_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value
