"""ME_PLOT input syntax, omitted defaults and photon-loop pitfalls."""

import pickle

import pytest

from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("MEPLOT")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_discovery_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen"
    assert not p.is_mpi()
    assert p.TASK.ITSEL() == p.TASK.NEPHOT() == 1
    assert p.TASK.INITIALSTATE() == "BAND"
    assert p.TASK.L() is None
    assert p.TASK.EPHOT() == 100.0
    assert p.TASK.CL() == "2P"
    assert p.TASK.SEFIN() == 0
    text = p.to_string(validate=True)
    assert "TASK\n\tMEPLOT" in text
    for name in ("ITSEL", "INITIALSTATE", "IPRINT", "L", "CL", "NONMAG", "SEFIN",
                 "EPHOT", "NEPHOT", "EPHOTMIN", "EPHOTMAX"):
        assert f"\t{name}=" not in p.TASK.to_string()
    assert "\tEMAX=" not in text
    assert p.ENERGY.EMAX() is None
    assert p.ENERGY.NE().tolist() == [300]
    assert p.ENERGY.GRID().tolist() == ["3"]
    assert p.ENERGY.ImE() == 0.01
    assert p.SITES.NL() is None
    path = tmp_path / "meplot.in"
    p.save_to_file(path, atoms=False)
    parsed = InputParameters.from_file(path)
    assert parsed.task_name == "meplot"
    assert parsed.to_dict() == p.to_dict()
    for copy in (p.copy(copy_values=True), pickle.loads(pickle.dumps(p))):
        assert copy.to_dict() == p.to_dict()
        copy.validate("save")


def test_iprint_inherits_control_and_uses_the_correct_key():
    p = parameters()
    p.CONTROL.PRINT = 2
    assert p.TASK.IPRINT() == 2
    assert "IPRINT=" not in p.TASK.to_string()
    p.TASK.IPRINT = 3
    text = p.to_string(validate=True)
    assert "\tIPRINT=3" in text
    assert p._definition.read_from_string(text).TASK.IPRINT() == 3


def test_explicit_default_nephot_is_not_written():
    p = parameters()
    p.TASK.NEPHOT = 1
    p.TASK.EPHOT = 150.0
    text = p.to_string(validate=True)
    assert "NEPHOT=" not in text
    assert "EPHOT=150.0" in text
    assert p._definition.read_from_string(text).TASK.EPHOT() == 150.0


def test_complete_photon_scan_and_relative_initial_energy():
    p = parameters()
    p.TASK.set({"NEPHOT": 5, "EPHOTMIN": 10.0, "EPHOTMAX": 100.0})
    p.ENERGY.set({"NE": 1, "EMINEV": -2.0, "EMAXEV": -2.0})
    text = p.to_string(validate=True)
    assert "NEPHOT=5" in text
    assert "EPHOTMIN=10.0" in text and "EPHOTMAX=100.0" in text
    assert "\tEPHOT=" not in text
    assert "EMINEV=-2.0" in text and "EMAXEV=-2.0" in text
    parsed = p._definition.read_from_string(text)
    assert parsed.to_dict() == p.to_dict()


@pytest.mark.parametrize("missing", ["EPHOTMIN", "EPHOTMAX"])
def test_photon_scan_requires_both_bounds(missing):
    p = parameters()
    p.ENERGY.NE = 1
    p.TASK.set({"NEPHOT": 5, "EPHOTMIN": 10.0, "EPHOTMAX": 100.0})
    p.TASK[missing].clear()
    with pytest.raises(DataValidityError, match="requires both"):
        p.validate("save")


@pytest.mark.parametrize("bounds", [{"EPHOTMIN": 10.0}, {"EPHOTMAX": 100.0},
                                  {"EPHOTMIN": 10.0, "EPHOTMAX": 100.0}])
def test_range_without_photon_scan_is_rejected(bounds):
    p = parameters()
    p.TASK.set(bounds)
    with pytest.raises(DataValidityError, match="require NEPHOT>1"):
        p.validate("save")


def test_photon_scan_cannot_silently_overwrite_initial_energy_results():
    p = parameters()
    p.TASK.set({"NEPHOT": 5, "EPHOTMIN": 10.0, "EPHOTMAX": 100.0})
    with pytest.raises(DataValidityError, match="overwrites"):
        p.validate("save")


def test_reversed_photon_range_is_rejected():
    p = parameters()
    p.ENERGY.NE = 1
    p.TASK.set({"NEPHOT": 5, "EPHOTMIN": 100.0, "EPHOTMAX": 10.0})
    with pytest.raises(DataValidityError, match="must be >="):
        p.validate("save")


def test_photon_scan_warns_about_ignored_sefin_and_fixed_energy():
    p = parameters()
    p.ENERGY.NE = 1
    p.TASK.set({"NEPHOT": 5, "EPHOTMIN": 10.0, "EPHOTMAX": 100.0,
                "EPHOT": 200.0, "SEFIN": 1})
    with pytest.warns(DataValidityWarning) as warnings:
        p.validate("save")
    assert any("disables TASK.SEFIN" in str(w.message) for w in warnings)
    assert any("EPHOT is ignored" in str(w.message) for w in warnings)


def test_core_uses_final_energy_grid_and_shell_without_photon_input():
    p = parameters()
    p.TASK.set({"INITIALSTATE": "CORE", "CL": "1S", "ITSEL": 2})
    p.ENERGY.set({"EMINEV": 0.0, "EMAXEV": 100.0})
    text = p.to_string(validate=True)
    assert "INITIALSTATE=CORE" in text and "CL=1S" in text and "ITSEL=2" in text
    assert p._definition.read_from_string(text).to_dict() == p.to_dict()
    p.TASK.EPHOT = 200.0
    with pytest.warns(DataValidityWarning, match="ignores TASK.EPHOT"):
        p.validate("save")


def test_core_shell_and_explicit_angular_cutoff_are_checked():
    p = parameters()
    p.TASK.INITIALSTATE = "CORE"
    p.TASK.CL = "1F"
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")
    p.TASK.CL = "4F"
    p.SITES.NL = 3
    with pytest.raises(DataValidityError, match="SITES.NL > l"):
        p.validate("save")
    p.SITES.NL = 4
    p.validate("save")


@pytest.mark.parametrize("name,value", [("ITSEL", 0), ("L", -1), ("NEPHOT", 0), ("SEFIN", 4)])
def test_invalid_parameter_values_are_rejected(name, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set(value)


def test_unverified_nonrelativistic_basis_warns():
    p = parameters()
    p.MODE.MODE = "NREL"
    with pytest.warns(DataValidityWarning, match="relativistic kappa"):
        p.validate("save")


def test_task_nonmag_is_distinct_from_control_nonmag():
    p = parameters()
    p.TASK.NONMAG = True
    text = p.to_string(validate=True)
    assert "NONMAG" in p.TASK.to_string()
    assert not p.CONTROL.NONMAG()
    assert p._definition.read_from_string(text).TASK.NONMAG()
