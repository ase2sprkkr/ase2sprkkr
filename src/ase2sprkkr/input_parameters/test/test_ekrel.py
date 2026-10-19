"""EKREL defaults, numbered paths and separation of search energies."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("EKREL")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen" and p.is_mpi()
    assert p.TASK.NK() == 1000 and p.TASK.NE() == 501
    assert p.TASK.EMIN() == 0.01 and p.TASK.EMAX() == 1
    assert p.TASK.ESLICE() is None and p.TASK.NKDIR() == 0
    text = p.to_string(validate=True)
    for key in ("NK", "NE", "EMIN", "EMAX", "ESLICE", "KPATH", "NKDIR"):
        assert f"\t{key}=" not in p.TASK.to_string()
    path = tmp_path / "ekrel.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_manual_path_round_trip():
    p = parameters()
    p.TASK.set({"NKDIR": 1, "KA1": [0.0, 0.0, 0.0], "KE1": [1.0, 0.0, 0.0]})
    text = p.to_string(validate=True)
    assert "KA1=" in text and "KE1=" in text
    assert p._definition.read_from_string(text).to_string() == text


def test_missing_manual_endpoint_rejected():
    p = parameters()
    p.TASK.NKDIR = 1
    with pytest.raises(DataValidityError, match="requires both"):
        p.validate("save")


def test_zero_length_segment_rejected():
    p = parameters()
    p.TASK.set({"NKDIR": 1, "KA1": [0.0, 0.0, 0.0], "KE1": [0.0, 0.0, 0.0]})
    with pytest.raises(DataValidityError, match="nonzero segment"):
        p.validate("save")


@pytest.mark.parametrize("key,value", [("NE", 1), ("NK", 1), ("NKDIR", 10)])
def test_invalid_counts(key, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[key].set(value)


def test_search_window_rejected():
    p = parameters()
    p.TASK.EMAX = 0.0
    with pytest.raises(DataValidityError, match="EMIN < TASK.EMAX"):
        p.validate("save")


def test_fermi_slice_round_trip():
    p = parameters()
    p.TASK['F-SURF'].set(True)
    p.TASK.ESLICE = 0.5
    text = p.to_string(validate=True)
    assert "F-SURF" in text and "ESLICE=0.5" in text
    assert p._definition.read_from_string(text).TASK.ESLICE() == 0.5
    p.TASK.ESLICE = 2.0
    with pytest.raises(DataValidityError, match="ESLICE"):
        p.validate("save")
