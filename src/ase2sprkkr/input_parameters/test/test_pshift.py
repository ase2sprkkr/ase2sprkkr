"""PSHIFT input round trips and warnings grounded in its source."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("PSHIFT")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.ENERGY.GRID().tolist() == ["1"]
    assert p.ENERGY.NE().tolist() == [300]
    assert p.ENERGY.ImE() == 0
    assert p.ENERGY.EMIN() == -0.2 and p.ENERGY.EMAX() == 1.0
    assert p.SITES.NL() is None
    text = p.to_string(validate=True)
    assert p.TASK.to_string().strip() == "TASK\n\tPSHIFT"
    assert "\tEMAX=1.0" in text and "\tNE={300}" in text
    path = tmp_path / "pshift.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL", "FREL"])
def test_relativistic_modes_round_trip(mode):
    p = parameters()
    p.MODE.MODE = mode
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).MODE.MODE() == p.MODE.MODE()


def test_relative_bounds_suppress_absolute_defaults():
    p = parameters()
    p.ENERGY.set({"EMINEV": -5.0, "EMAXEV": 10.0, "NE": 80})
    text = p.to_string(validate=True)
    assert "EMINEV=-5.0" in text and "EMAXEV=10.0" in text
    assert "\tEMAX=" not in text and "\tEMIN=" not in text
    assert p._definition.read_from_string(text).to_string() == text


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_unsupported_correlations(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_imaginary_energy_is_overwritten():
    p = parameters()
    p.ENERGY.ImE = 0.01
    with pytest.warns(DataValidityWarning, match="resets ENERGY.ImE"):
        p.validate("save")


def test_split_and_second_path_warn():
    p = parameters()
    p.ENERGY.set({"GRID": [5, 3], "NE": [20, 30], "SPLITSS": True})
    with pytest.warns(DataValidityWarning) as records:
        p.validate("save")
    messages = [str(w.message) for w in records]
    assert any("disables ENERGY.SPLITSS" in m for m in messages)
    assert any("only the first ENERGY path" in m for m in messages)
    assert any("unwraps phases" in m for m in messages)


def test_relativistic_soc_scaling_warns():
    p = parameters()
    p.MODE.SOC = 0.5
    with pytest.warns(DataValidityWarning, match="without SOCTL"):
        p.validate("save")


def test_nonuniform_nonrelativistic_cutoff_warns():
    p = parameters()
    p.MODE.MODE = "SP-SREL"
    p.SITES.NL = [2, 3]
    with pytest.warns(DataValidityWarning, match="uninitialized"):
        p.validate("save")


def test_high_relativistic_cutoff_warns_about_legends():
    p = parameters()
    p.SITES.NL = 6
    with pytest.warns(DataValidityWarning, match="invalid legends"):
        p.validate("save")


@pytest.mark.parametrize("ne", [0, -1])
def test_nonpositive_energy_count(ne):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = ne
