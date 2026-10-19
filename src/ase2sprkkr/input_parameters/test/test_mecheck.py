"""MECHECK syntax and source-specific diagnostic restrictions."""
import pickle
import numpy as np
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("MECHECK")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen" and not p.is_mpi()
    assert p.TASK.ERYDI().tolist() == [0.5, 0.01]
    assert p.TASK.ERYDF().tolist() == [0.8, 0.01]
    assert p.SITES.NL() is None and p.ENERGY.EMAX() is None
    text = p.to_string(validate=True)
    assert "TASK\n\tMECHECK" in text
    for key in ("ITSEL", "PRINT", "INITIALSTATE", "CL", "ERYDI", "ERYDF", "NMEFORM"):
        assert f"\t{key}=" not in p.TASK.to_string()
    path = tmp_path / "mecheck.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_dict() == p.to_dict()
    assert pickle.loads(pickle.dumps(p)).to_string() == p.to_string()
    assert p.copy(copy_values=True).to_string() == p.to_string()


def test_complex_energies_and_print():
    p = parameters()
    p.CONTROL.PRINT = 2
    assert p.TASK.PRINT() == 2
    p.TASK.set({"ERYDI": [-0.1, 0.02], "ERYDF": [0.9, 0.03], "PRINT": 3, "ITSEL": 2})
    text = p.to_string(validate=True)
    assert "\tPRINT=3" in p.TASK.to_string() and "IPRINT" not in text
    parsed = p._definition.read_from_string(text)
    assert parsed.to_string() == p.to_string()
    np.testing.assert_array_equal(parsed.TASK.ERYDI(), p.TASK.ERYDI())
    np.testing.assert_array_equal(parsed.TASK.ERYDF(), p.TASK.ERYDF())


@pytest.mark.parametrize("key", ["ERYDI", "ERYDF"])
def test_complex_energy_requires_two_components(key):
    with pytest.raises(DataValidityError):
        parameters().TASK[key].set([0.5])


@pytest.mark.parametrize("key", ["ERYDI", "ERYDF"])
def test_zero_energy_warns(key):
    p = parameters()
    p.TASK[key].set([0.0, 0.0])
    with pytest.warns(DataValidityWarning, match="magnitude"):
        p.validate("save")


@pytest.mark.parametrize("value", [1, 2, 3])
def test_nmeform_is_overwritten(value):
    p = parameters()
    p.TASK.NMEFORM = value
    with pytest.warns(DataValidityWarning, match="overwrites"):
        p.validate("save")


def test_core_scan_and_ignored_energies():
    p = parameters()
    p.TASK.set({"INITIALSTATE": "CORE", "CL": "1S"})
    p.ENERGY.set({"NE": 20, "EMINEV": 0.0, "EMAXEV": 100.0})
    text = p.to_string(validate=True)
    parsed = p._definition.read_from_string(text)
    assert parsed.to_string() == p.to_string()
    np.testing.assert_array_equal(parsed.TASK.ERYDI(), p.TASK.ERYDI())
    np.testing.assert_array_equal(parsed.TASK.ERYDF(), p.TASK.ERYDF())
    p.TASK.ERYDF = [0.9, 0.02]
    with pytest.warns(DataValidityWarning, match="ERYDF is ignored"):
        p.validate("save")


def test_core_shell_validation():
    p = parameters()
    p.TASK.set({"INITIALSTATE": "CORE", "CL": "1P"})
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")
    p.TASK.CL = "4F"
    p.SITES.NL = 3
    with pytest.raises(DataValidityError, match="SITES.NL"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_relativistic_requirement(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="IREL=3"):
        p.validate("save")


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_unsupported_correlations(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_band_warns_about_energy_grid():
    p = parameters()
    p.ENERGY.NE = 10
    with pytest.warns(DataValidityWarning, match="ENERGY does not select"):
        p.validate("save")


def test_task_nonmag_is_distinct_from_control():
    p = parameters()
    p.TASK.NONMAG = True
    text = p.to_string(validate=True)
    assert "NONMAG" in p.TASK.to_string() and not p.CONTROL.NONMAG()
    assert p._definition.read_from_string(text).TASK.NONMAG()
