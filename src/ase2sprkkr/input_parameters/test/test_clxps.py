"""Core-state allocation, native geometry names and omitted defaults."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("CLXPS")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.executable == "kkrgen" and not p.task_definition.mpi
    assert p.TASK.CL() == "2P" and p.TASK.EPHOT() == 1253.6
    assert p.TASK.ImE() == 0.01 and p.ENERGY.NE().tolist() == [18]
    text = p.to_string(validate=True)
    for key in ("ME", "IT", "CL", "EPHOT", "ImE"):
        assert f"\t{key}=" not in p.TASK.to_string()
    path = tmp_path / "clxps.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_geometry_zero_round_trip():
    p = parameters()
    p.TASK.set({"ANGRES": True, "ANGQ0": [0.0, 0.0], "ANGE0": [90.0, 0.0], "ANGK0": [45.0, 0.0]})
    text = p.to_string(validate=True)
    assert "ANGQ0=" in text and "ANGK0=" in text
    assert p._definition.read_from_string(text).to_string() == text


def test_incomplete_geometry_rejected():
    p = parameters()
    p.TASK.ANGQ1 = [0.0, 0.0]
    with pytest.raises(DataValidityError, match="requires ANGQ1"):
        p.validate("save")


def test_gap_in_geometries_rejected():
    p = parameters()
    p.TASK.set({"ANGQ2": [0.0, 0.0], "ANGE2": [90.0, 0.0], "ANGK2": [45.0, 0.0]})
    with pytest.raises(DataValidityError, match="consecutive"):
        p.validate("save")


def test_angle_resolution_needs_geometry():
    p = parameters()
    p.TASK.ANGRES = True
    with pytest.raises(DataValidityError, match="complete geometry"):
        p.validate("save")


def test_core_allocation_and_shell():
    p = parameters()
    p.TASK.CL = "4F"
    p.ENERGY.NE = 6
    with pytest.raises(DataValidityError, match="4.l.2"):
        p.validate("save")
    p.ENERGY.NE = 18
    p.validate("save")
    p.TASK.CL = "1P"
    with pytest.raises(DataValidityError, match="l < n"):
        p.validate("save")


def test_zero_imaginary_energy_warns():
    p = parameters()
    p.TASK.ImE = 0.0
    with pytest.warns(DataValidityWarning, match="replaces near-zero"):
        p.validate("save")


@pytest.mark.parametrize("mode", ["NREL", "SREL", "SP-SREL"])
def test_relativistic_requirement(mode):
    p = parameters()
    p.MODE.MODE = mode
    with pytest.raises(DataValidityError, match="IREL=3"):
        p.validate("save")


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_correlations_rejected(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


def test_unimplemented_matrix_element_form_rejected():
    with pytest.raises(DataValidityError):
        parameters().TASK.ME = "GRV"
