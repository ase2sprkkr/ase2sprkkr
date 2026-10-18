import io
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters():
    p = InputParameters.create("CLPES")
    p.CONTROL.POTFIL.set("Fe.pot")
    return p


def test_roundtrip():
    p = parameters()
    q = p._definition.read_from_file(io.StringIO(p.to_string(validate=True)))
    assert q.TASK.TASK() == "CLPES"
    assert q.SPEC_CL.CL() == "2P"
    assert q.SPEC_CL.IT() == 1
    assert not q.task_definition.mpi


@pytest.mark.parametrize("option,value", [("SPEC_CL.CL", "1P"), ("SPEC_CL.IT", 0), ("SPEC_PH.NPOL", 0), ("SPEC_PH.NPOL", 3), ("SPEC_EL.NT", 2), ("ENERGY.EWORK_EV", -1.)])
def test_invalid(option, value):
    p = parameters()
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")


def test_single_core_state_bounds():
    p = parameters()
    p.ENERGY.NE.set([1])
    p.SPEC_CL.IE_ICST.set(7)
    with pytest.raises((DataValidityError, DataValidityErrors), match="state count"):
        p.validate("save")
