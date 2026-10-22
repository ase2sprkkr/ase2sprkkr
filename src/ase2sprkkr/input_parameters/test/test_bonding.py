import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("BONDING")
    p.CONTROL.POTFIL = "molecule.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert not p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    assert p.TAU.MOL() and p.TAU.BZINT() == "CLUSTER" and p.TAU.NSHLCLU() == 0
    text = p.to_string(validate=True)
    assert "EMAX=1.0" in text and "MOL" in text and "NSHLCLU=0" in text
    path = tmp_path / "bonding.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


@pytest.mark.parametrize("section,name,value", [
    ("TAU", "BZINT", "POINTS"), ("TAU", "MOL", False),
    ("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"), ("MODE", "MODE", "SP-SREL"),
    ("ENERGY", "SPLITSS", True), ("ENERGY", "GRID", [3, 5]), ("ENERGY", "NE", [300, 30]),
])
def test_restrictions(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("value", [0, 1, -1])
def test_invalid_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value
