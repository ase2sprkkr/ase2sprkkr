import io
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def test_tdobs_roundtrip():
    p = InputParameters.create("TDOBS")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.TAU.TAUDIR.set("old.tau.dir")
    p.TAU.TAUPATH.set("cache")
    q = p._definition.read_from_file(io.StringIO(p.to_string(validate=True)))
    assert q.TAU.TAUDIR() == "old.tau.dir"
    assert q.TAU.TAUPATH() == "cache"
    assert q.TASK.TASK() == "TDOBS"
    assert q.CONTROL.ADSI() == "TDOBS"
    assert not q.task_definition.mpi
    assert q.task_definition.executable == "kkrgen"


@pytest.mark.parametrize("option,value", [("ENERGY.NE", [0]), ("ENERGY.NE", [1]), ("ENERGY.ImE", -1.), ("ENERGY.EMIN", float("nan"))])
def test_tdobs_invalid_mesh(option, value):
    p = InputParameters.create("TDOBS")
    p.CONTROL.POTFIL.set("Fe.pot")
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")
