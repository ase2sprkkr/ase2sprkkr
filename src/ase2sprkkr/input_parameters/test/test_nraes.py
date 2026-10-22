import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("NRAES")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert p.MODE.MODE() == "SREL" and not p.task_definition.mpi
    assert p.TASK.CL() == "2P"
    text = p.to_string(validate=True)
    for name in ("IT", "CL"):
        assert f"{name}=" not in p.TASK.to_string()
    path = tmp_path / "nraes.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_explicit_options_round_trip():
    p = parameters()
    p.TASK.set({"IT": 2, "CL": "3D", "CALCLS": True, "OUTGNU": True})
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).to_string() == text


def test_invalid_shell():
    p = parameters()
    p.TASK.CL = "1P"
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_correlations_rejected(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_mode_override_warning():
    p = parameters()
    p.MODE.MODE = "SP-SREL"
    with pytest.warns(DataValidityWarning, match="forces IREL"):
        p.validate("save")


def test_capacity_rejected():
    p = parameters()
    p.ENERGY.NE = 1
    with pytest.raises(DataValidityError, match="two DOS energies"):
        p.validate("save")


def test_shell_label_table_limit():
    with pytest.raises(DataValidityError):
        parameters().TASK.CL = "6S"
