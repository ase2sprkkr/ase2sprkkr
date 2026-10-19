"""SOCPAR discovery, real-energy input and overwritten settings."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("SOCPAR")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen" and not p.is_mpi()
    assert p.ENERGY.GRID().tolist() == ["1"]
    assert p.ENERGY.NE().tolist() == [300]
    assert p.ENERGY.ImE() == 0
    assert p.ENERGY.EMAX() is None and p.SITES.NL() is None
    text = p.to_string(validate=True)
    assert p.TASK.to_string().strip() == "TASK\n\tSOCPAR"
    assert "\tGRID={1}" in text and "\tNE={300}" in text
    assert "\tEMAX=" not in text
    path = tmp_path / "socpar.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


def test_relative_energy_window():
    p = parameters()
    p.ENERGY.set({"EMINEV": -5.0, "EMAXEV": 10.0, "NE": 80})
    text = p.to_string(validate=True)
    assert "EMINEV=-5.0" in text and "EMAXEV=10.0" in text
    assert p._definition.read_from_string(text).to_string() == text


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_unsupported_correlations(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_unverified_relativistic_modes_warn(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.warns(DataValidityWarning, match="relativistic wave-function"):
        p.validate("save")


def test_soc_scaling_is_overwritten():
    p = parameters()
    p.MODE.SOC = 0.5
    with pytest.warns(DataValidityWarning, match="SOCTL=1e-3"):
        p.validate("save")


def test_imaginary_energy_is_overwritten():
    p = parameters()
    p.ENERGY.ImE = 0.01
    with pytest.warns(DataValidityWarning, match="resets ENERGY.ImE"):
        p.validate("save")


def test_second_path_and_split_are_ignored():
    p = parameters()
    p.ENERGY.set({"GRID": [5, 3], "NE": [20, 30], "SPLITSS": True})
    with pytest.warns(DataValidityWarning) as records:
        p.validate("save")
    messages = [str(w.message) for w in records]
    assert any("disables ENERGY.SPLITSS" in m for m in messages)
    assert any("only the first ENERGY path" in m for m in messages)
    assert any("projected onto the real axis" in m for m in messages)


@pytest.mark.parametrize("nl", [1, 52])
def test_invalid_orbital_cutoff(nl):
    p = parameters()
    p.SITES.NL = nl
    with pytest.raises(DataValidityError, match="SOCPAR"):
        p.validate("save")


@pytest.mark.parametrize("nl", [2, 3, 4])
def test_valid_orbital_cutoff(nl):
    p = parameters()
    p.SITES.NL = nl
    p.validate("save")


@pytest.mark.parametrize("ne", [0, -1])
def test_nonpositive_energy_count(ne):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = ne
