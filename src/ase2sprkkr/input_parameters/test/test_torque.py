import pickle
import numpy as np
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("TORQUE")
    p.CONTROL.POTFIL = "rotated.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert p.is_mpi() and p._definition.executable == "kkrgen"
    assert p.ENERGY.GRID().tolist() == ["8"] and p.ENERGY.NE().tolist() == [36]
    assert p.TASK.THETAQ().tolist() == [90.0]
    text = p.to_string(validate=True)
    path = tmp_path / "torque.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_site_angles(tmp_path):
    p = parameters()
    p.TASK.THETAQ = [90., 45.]
    p.TASK.PHIQ = [0., 90.]
    path = tmp_path / "sites.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.THETAQ().tolist() == [90., 45.]
    assert q.TASK.PHIQ().tolist() == [0., 90.]


@pytest.mark.parametrize("name", ["THETAQ", "PHIQ"])
def test_empty_angles(name):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set([])


@pytest.mark.parametrize("name,value", [("THETAQ", np.nan), ("PHIQ", np.inf)])
def test_nonfinite_angles(name, value):
    p = parameters()
    p.TASK[name].set([value])
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("SPLITSS", True), ("GRID", [8, 3]), ("NE", [36, 10])])
def test_split_rejected(name, value):
    p = parameters()
    p.ENERGY[name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("value", [0, -1])
def test_invalid_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value
