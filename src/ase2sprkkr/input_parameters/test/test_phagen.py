import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("PHAGEN")
    p.CONTROL.POTFIL = "Cu.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert not p.task_definition.mpi and p.task_definition.executable == "kkrscf"
    assert p.SCF.VXC() == "VWN"
    assert p.ENERGY.NE().tolist() == [32]
    text = p.to_string(validate=True)
    assert "PHAGEN_CALCTYPE" not in p.TASK and "PHAGEN" in text
    path = tmp_path / "phagen.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


@pytest.mark.parametrize("section,name,value", [
    ("SCF", "FULLPOT", True), ("SCF", "VXC", "PBE"), ("SCF", "VXC", "PW92"),
    ("MODE", "OP", "DMFT"), ("MODE", "OP", "LDA+U"),
    ("ENERGY", "SPLITSS", True), ("ENERGY", "GRID", [5, 3]), ("ENERGY", "NE", [32, 10]),
])
def test_unsupported_modes(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("value", [0, -1])
def test_invalid_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value


@pytest.mark.parametrize("mode", ["NREL", "SREL"])
def test_no_invented_relativity_ban(mode):
    p = parameters()
    p.MODE.MODE = mode
    p.validate("save")
