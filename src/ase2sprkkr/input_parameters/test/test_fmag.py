import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("FMAG")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    assert p.ENERGY.GRID().tolist() == ["5"] and p.ENERGY.NE().tolist() == [30]
    text = p.to_string(validate=True)
    assert "EMAX=" not in text
    assert "FMAG" in p.TASK.to_string()
    path = tmp_path / "fmag.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


@pytest.mark.parametrize("name,value", [("SPLITSS", True), ("GRID", [5, 3]), ("NE", [30, 10])])
def test_split_contour_rejected(name, value):
    p = parameters()
    p.ENERGY[name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


def test_other_contour_warns():
    p = parameters()
    p.ENERGY.GRID = 3
    with pytest.warns(DataValidityWarning, match="GRID=5"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL"])
def test_nonmagnetic_normalization_warning(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.warns(DataValidityWarning, match="divide by zero"):
        p.validate("save")


def test_empty_mesh_rejected():
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = 0
