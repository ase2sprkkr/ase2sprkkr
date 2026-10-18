"""Regression tests for the EF workflow and native Stoner constraints."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("STONER")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert not p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    assert p.ENERGY.GRID().tolist() == ["0"]
    assert p.ENERGY.NE().tolist() == [1]
    assert p.ENERGY.ImE() == 0
    text = p.to_string(validate=True)
    assert "EMIN=" not in text and "EMAX=" not in text
    path = tmp_path / "stoner.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


@pytest.mark.parametrize("name,value", [("GRID", 3), ("NE", 2), ("ImE", .01),
                                         ("SPLITSS", True), ("EMIN", -.2), ("EMAXEV", 1.)])
def test_energy_modifications_rejected(name, value):
    p = parameters()
    p.ENERGY[name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


def test_relative_zero_is_fermi_energy():
    p = parameters()
    p.ENERGY.EMINEV = 0.0
    p.ENERGY.EMAXEV = 0.0
    p.to_string(validate=True)


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_nonrelativistic_wave_function_format_rejected(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="WAVFUN_READ_REL"):
        p.validate("save")


def test_missing_d_channel_rejected():
    p = parameters()
    p.SITES.NL = [2]
    with pytest.raises(DataValidityError, match="d channel"):
        p.validate("save")
    p.SITES.NL = [3]
    p.validate("save")
