"""PLOTPRS round trips and initialization requirements."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("PLOTPRS")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen" and not p.is_mpi()
    assert p.ENERGY.NE().tolist() == [1]
    assert p.ENERGY.GRID().tolist() == ["1"]
    assert p.ENERGY.EMAX() is None and p.SITES.NL() is None
    text = p.to_string(validate=True)
    assert p.TASK.to_string().strip() == "TASK\n\tPLOTPRS"
    assert "\tNE={1}" in text and "\tEMAX=" not in text
    path = tmp_path / "plotprs.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_unsupported_correlations(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL", "FREL"])
def test_no_invented_relativistic_restriction(mode):
    p = parameters()
    p.MODE.MODE = mode
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).MODE.MODE() == p.MODE.MODE()


def test_hyperfine_flag_is_ignored():
    p = parameters()
    p.TASK.HFF = True
    with pytest.warns(DataValidityWarning, match="hyperfine"):
        p.validate("save")


def test_energy_count_does_not_define_plot_samples():
    p = parameters()
    p.ENERGY.NE = 30
    with pytest.warns(DataValidityWarning, match="radial data"):
        p.validate("save")


@pytest.mark.parametrize("ne", [0, -1])
def test_initialization_needs_positive_energy_count(ne):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = ne


def test_explicit_shape_function_filename_round_trip():
    p = parameters()
    p.CONTROL.SFNFIL = "custom.sfn"
    text = p.to_string(validate=True)
    assert "SFNFIL=custom.sfn" in text
    assert p._definition.read_from_string(text).CONTROL.SFNFIL() == "custom.sfn"
