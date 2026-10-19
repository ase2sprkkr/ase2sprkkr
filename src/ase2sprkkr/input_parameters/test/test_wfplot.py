"""WFPLOT reader syntax, omitted defaults and source-level restrictions."""

import pickle
import re
from contextlib import nullcontext

import pytest

from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    result = InputParameters.create("WFPLOT")
    result.CONTROL.POTFIL = "Fe.pot"
    return result


def test_discovery_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen"
    assert not p.is_mpi()
    text = p.to_string(validate=True)
    assert "TASK\n\tWFPLOT" in text
    assert "\tITSEL=1" in text
    assert "\tSTATE=" not in text
    assert "\tCL=" not in text
    assert "\tMJ=" not in text
    assert "\tL=" not in text
    assert p.TASK.STATE() == "BAND"
    assert p.TASK.CL() == "2P"
    assert p.TASK.MJ() == "1/2"
    assert p.TASK.L() is None
    assert p.ENERGY.NE().tolist() == [1]
    assert re.search(r"\tGRID=\{?3\}?\n", text)
    assert re.search(r"\tNE=\{?1\}?\n", text)
    assert p.ENERGY.EMIN() == p.ENERGY["EMAX"]() == -0.2
    assert p.ENERGY.ImE() == 0.01
    path = tmp_path / "wfplot.in"
    p.save_to_file(path, atoms=False)
    parsed = InputParameters.from_file(path)
    assert parsed.task_name == "wfplot"
    assert parsed.to_dict() == p.to_dict()
    for copied in (p.copy(copy_values=True), pickle.loads(pickle.dumps(p))):
        assert copied.to_dict() == p.to_dict()
        copied.validate("save")


def test_print_inherits_control_and_can_be_overridden():
    p = parameters()
    p.CONTROL.PRINT = 2
    assert p.TASK.PRINT() == 2
    assert "PRINT=" not in p.TASK.to_string()
    p.TASK.PRINT = 3
    assert "PRINT=3" in p.TASK.to_string()
    parsed = p._definition.read_from_string(p.to_string(validate=True))
    assert parsed.CONTROL.PRINT() == 2 and parsed.TASK.PRINT() == 3


@pytest.mark.parametrize("value", ["1/2", "+3/2", "-5/2", "+61/2", "-61/2"])
def test_mj_uses_fortran_fraction_syntax(value):
    p = parameters()
    p.TASK.MJ = value
    expected = pytest.warns(DataValidityWarning, match="reduces TASK.MJ") if "61" in value else nullcontext()
    with expected:
        parsed = p._definition.read_from_string(p.to_string(validate=False))
    assert parsed.TASK.MJ() == value


@pytest.mark.parametrize("value", [0.5, "0.5", "2/2", "63/2", "1/3"])
def test_invalid_mj_is_rejected(value):
    with pytest.raises(DataValidityError):
        parameters().TASK.MJ = value


def test_core_selection_and_source_mj_clamping():
    p = parameters()
    p.TASK.set({"STATE": "CORE", "CL": "1S", "MJ": "3/2", "ITSEL": 2})
    with pytest.warns(DataValidityWarning, match="reduces TASK.MJ"):
        text = p.to_string(validate=True)
    assert "STATE=CORE" in text and "CL=1S" in text and "ITSEL=2" in text
    p.TASK.MJ = "-1/2"
    p.validate("save")
    p.TASK.CL = "1F"
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")


def test_band_l_tracks_uniform_cutoff_without_serializing_a_fake_default():
    p = parameters()
    p.SITES.NL = 4
    assert p.TASK.L() is None
    p.TASK.L = 4
    with pytest.raises(DataValidityError, match="smaller than SITES.NL"):
        p.validate("save")
    p.TASK.L = 3
    p.validate("save")


def test_relative_single_point_energy_writes_both_bounds():
    p = parameters()
    p.ENERGY.EMINEV = 2.0
    text = p.to_string(validate=True)
    assert "\tEMINEV=2.0" in text and "\tEMAXEV=2.0" in text
    assert "\tEMIN=" not in text and "\tEMAX=" not in text


@pytest.mark.parametrize("method", ["DMFT", "LDA+U"])
def test_unsupported_many_body_methods_are_rejected(method):
    p = parameters()
    p.MODE.OP = method
    with pytest.raises(DataValidityError, match="does not support MODE.OP"):
        p.validate("save")


def test_ignored_mode_and_extra_energy_points_warn():
    p = parameters()
    p.MODE.MODE = "NREL"
    with pytest.warns(DataValidityWarning, match="forces fully relativistic"):
        p.validate("save")
    p.MODE.MODE = "FREL"
    p.ENERGY.NE = 3
    with pytest.warns(DataValidityWarning, match="only the first energy point"):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("ITSEL", 0), ("L", -1)])
def test_invalid_indices_are_rejected(name, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set(value)
