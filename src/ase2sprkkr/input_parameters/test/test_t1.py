"""T1 discovery and input safeguards against native initialization defects."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("T1")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.CONTROL.PRINT() == 1 and p.CONTROL.RUNELOOP()
    assert p.ENERGY.ImE() == 0 and p.ENERGY.NE().tolist() == [1]
    assert p.SITES.NL() is None
    text = p.to_string(validate=True)
    assert "\tPRINT=1" in text and "\tRUNELOOP" in text
    assert "\tEMIN=" not in text and "\tEMAX=" not in text
    path = tmp_path / "t1.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


@pytest.mark.parametrize("mode", ["NREL", "SREL"])
def test_disabled_nonrelativistic_implementation(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="not available"):
        p.validate("save")


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_correlations_rejected(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_zero_print_rejected():
    p = parameters()
    p.CONTROL.PRINT = 0
    with pytest.raises(DataValidityError, match="NN/MM"):
        p.validate("save")


def test_disabled_loop_rejected():
    p = parameters()
    p.CONTROL.RUNELOOP = False
    with pytest.raises(DataValidityError, match="initialize TAUT"):
        p.validate("save")


def test_energy_grid_is_overwritten():
    p = parameters()
    p.ENERGY.set({"GRID": 1, "NE": 20})
    with pytest.warns(DataValidityWarning, match="overwritten"):
        p.validate("save")


def test_relativistic_label_capacity():
    p = parameters()
    p.SITES.NL = 6
    with pytest.raises(DataValidityError, match="10 kappa labels"):
        p.validate("save")


def test_imaginary_fermi_energy_round_trip():
    p = parameters()
    p.ENERGY.ImE = 0.01
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).ENERGY.ImE() == 0.01
