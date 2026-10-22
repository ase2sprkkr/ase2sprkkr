import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("AES")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.TASK.NEME() == 10 and p.TASK.CL() == "2P"
    text = p.to_string(validate=True)
    for name in ("NEME", "IT", "CL"):
        assert f"{name}=" not in p.TASK.to_string()
    assert "EMAX=" not in text
    path = tmp_path / "aes.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_explicit_options_round_trip():
    p = parameters()
    p.TASK.set({"IT": 2, "CL": "3D", "NEME": 12, "PRINTME": True, "PRINTAME": True})
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).to_string() == text


@pytest.mark.parametrize("points", [0, 1, 2])
def test_interpolation_mesh_minimum(points):
    with pytest.raises(DataValidityError):
        parameters().TASK.NEME = points


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_relativistic_requirement(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="IREL=3"):
        p.validate("save")


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_correlations_rejected(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_invalid_shell():
    p = parameters()
    p.TASK.CL = "1P"
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")
    with pytest.raises(DataValidityError):
        p.TASK.CL = "6S"


def test_capacity_minimum():
    p = parameters()
    p.ENERGY.NE = 2
    with pytest.raises(DataValidityError, match="three TAU energies"):
        p.validate("save")


def test_core_cache_name_table_limit():
    p = parameters()
    p.TASK.CL = "4F"
    with pytest.raises(DataValidityError, match="STRME"):
        p.validate("save")
