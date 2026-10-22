import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("VBPES")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.TASK.EPHOT() == 1253.6 and p.TASK.MODE() == 2
    assert p.TASK.ImEFIN() is None and p.TASK.NSPIN() is None
    text = p.to_string(validate=True)
    assert "ImEFIN=" not in text and "EMAX=" not in text
    assert "EPHOT=" not in p.TASK.to_string()
    path = tmp_path / "vbpes.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_custom_geometry_round_trip():
    p = parameters()
    p.TASK.set({"GEOMETRY": 3, "TETPH": 45.0, "PHIPH": 90.0, "ImEFIN": .02})
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).to_string() == text


def test_custom_vector_geometry():
    p = parameters()
    p.TASK.set({"GEOMETRY": 4, "QVECPH": [1.0, 0.0, 1.0]})
    p.validate("save")


@pytest.mark.parametrize("name", ["QVEC", "QVECPH", "NSPIN"])
def test_zero_vector_rejected(name):
    p = parameters()
    p.TASK[name].set([0.0, 0.0, 0.0])
    with pytest.raises(DataValidityError, match="nonzero"):
        p.validate("save")


def test_custom_geometry_needs_direction():
    p = parameters()
    p.TASK.GEOMETRY = 3
    with pytest.raises(DataValidityError, match="Custom geometry"):
        p.validate("save")
    p.TASK.TETPH = 45.0
    with pytest.raises(DataValidityError, match="requires PHIPH"):
        p.validate("save")


def test_final_state_flag_precedence():
    p = parameters()
    p.TASK.set({"USETFIN": True, "USETAUNN": True})
    with pytest.warns(DataValidityWarning, match="disables USETAUNN"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_relativistic_requirement(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="IREL=3"):
        p.validate("save")


def test_grv_rejected():
    with pytest.raises(DataValidityError):
        parameters().TASK.ME = "GRV"
