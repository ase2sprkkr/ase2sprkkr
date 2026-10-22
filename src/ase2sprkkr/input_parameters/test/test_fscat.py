"""FSCAT input round trips and meaningful source restrictions."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("FSCAT")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.ENERGY.GRID().tolist() == ["1"]
    assert p.ENERGY.NE().tolist() == [300]
    assert p.ENERGY.ImE() == 0
    assert p.ENERGY.EMAX() is None and p.SITES.NL() is None
    text = p.to_string(validate=True)
    assert p.TASK.to_string().strip() == "TASK\n\tFSCAT"
    assert "\tEMAX=" not in text
    path = tmp_path / "fscat.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


def test_relative_window():
    p = parameters()
    p.ENERGY.set({"EMINEV": 0.0, "EMAXEV": 100.0, "NE": 80})
    text = p.to_string(validate=True)
    assert "EMINEV=0.0" in text and "EMAXEV=100.0" in text
    assert p._definition.read_from_string(text).to_string() == text


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_unsupported_correlations(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_unverified_basis_warns(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.warns(DataValidityWarning, match="kappa/spin"):
        p.validate("save")


def test_overwritten_imaginary_energy_warns():
    p = parameters()
    p.ENERGY.ImE = 0.01
    with pytest.warns(DataValidityWarning, match="resets ENERGY.ImE"):
        p.validate("save")


def test_second_split_path_warns():
    p = parameters()
    p.ENERGY.set({"GRID": [5, 3], "NE": [20, 30], "SPLITSS": True})
    with pytest.warns(DataValidityWarning) as records:
        p.validate("save")
    assert any("disables ENERGY.SPLITSS" in str(w.message) for w in records)
    assert any("only the first ENERGY path" in str(w.message) for w in records)


def test_s_only_cutoff_rejected():
    p = parameters()
    p.SITES.NL = 1
    with pytest.raises(DataValidityError, match="CALC_CHPLM"):
        p.validate("save")


def test_nonuniform_cutoffs_warn():
    p = parameters()
    p.SITES.NL = [2, 3]
    with pytest.warns(DataValidityWarning, match="global maximum"):
        p.validate("save")


@pytest.mark.parametrize("ne", [0, -1])
def test_invalid_energy_count(ne):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = ne
