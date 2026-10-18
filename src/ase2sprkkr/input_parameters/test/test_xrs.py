import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("XRS")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert not p.is_mpi() and p._definition.executable == "kkrgen"
    assert p.TASK.NE3() == 40 and p.TASK.EMAX3() == 8.0 and p.TASK.TSELECT() == 0
    text = p.to_string(validate=True)
    assert "EMIN=" not in text and "NE3=" not in text and "EMAX=4.0" in text
    assert "OUTPUT" not in p.TASK
    path = tmp_path / "xrs.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_resonance_options(tmp_path):
    p = parameters()
    p.TASK.set({"NE3": 60, "EMAX3": 10.0, "TSELECT": 1, "TAUCORE": [.3, .5], "POLOBS": "C"})
    path = tmp_path / "options.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.NE3() == 60 and q.TASK.TSELECT() == 1
    assert q.TASK.TAUCORE().tolist() == [.3, .5]


@pytest.mark.parametrize("section,name,value", [
    ("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"), ("MODE", "MODE", "SP-SREL"),
    ("ENERGY", "GRID", [6, 3]), ("ENERGY", "NE", [180, 30]), ("ENERGY", "SPLITSS", True),
    ("TASK", "CL", "1P"), ("TASK", "TAUCORE", [-.1]),
    ("TASK", "EMAX3", 4.0), ("TASK", "EMAX3", float("inf")),
])
def test_invalid_combinations(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("IT", 0), ("CL", "6S"), ("NE3", 0), ("TSELECT", 3)])
def test_invalid_values(name, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set(value)
