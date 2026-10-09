"""SCF native algorithm and input-domain regressions."""
import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("SCF")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


@pytest.mark.parametrize("algorithm", ["BROYDEN2", "ANDERSON", "TCHEBY"])
def test_algorithms_roundtrip(algorithm):
    p = parameters()
    p.SCF.ALG.set(algorithm)
    q = p._definition.read_from_file(io.StringIO(p.to_string()))
    assert q.SCF.ALG() == algorithm
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    assert p._definition.executable == "kkrscf" and p._definition.mpi is True


def test_defaults_preserved():
    p = parameters()
    assert p.SCF.NITER() == 200 and p.SCF.MIX() == 0.2
    assert p.SCF.TOL() == 1e-5 and p.SCF.ITDEPT() == 40
    assert p.SCF.ALG() == "BROYDEN2"
    assert p.ENERGY.GRID().tolist() == ["5"]


@pytest.mark.parametrize("option,value", [
    ("SCF.ALG", "PULAY"), ("SCF.NITER", 0),
    ("SCF.ISTBRY", 0), ("SCF.ITDEPT", 2), ("SCF.ITDEPT", 41),
    ("SCF.MIX", float("nan")), ("SCF.TOL", 0.),
    ("SCF.EFGUESS", float("inf")), ("ENERGY.NE", [0]),
])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")


def test_single_point_rectangular_mesh():
    p = parameters()
    p.ENERGY.GRID.set([3])
    p.ENERGY.NE.set([1])
    with pytest.raises((DataValidityError, DataValidityErrors), match="NE>=2"):
        p.validate("save")


def test_tcheby_fullpot():
    p = parameters()
    p.SCF.ALG.set("TCHEBY")
    p.SCF.FULLPOT.set(True)
    with pytest.raises((DataValidityError, DataValidityErrors), match="TCHEBY"):
        p.validate("save")
