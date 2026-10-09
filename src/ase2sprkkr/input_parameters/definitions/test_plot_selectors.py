import io
import pickle
import pytest
from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError, DataValidityErrors


@pytest.mark.parametrize("task", ["PLOTPOT", "PLOTRHO", "PLOTSFN"])
def test_roundtrip_and_pickle(task):
    p = InputParameters.create(task)
    p.CONTROL.POTFIL.set("Fe.pot")
    q = p._definition.read_from_file(io.StringIO(p.to_string(validate=True)))
    assert q.TASK.TASK() == task
    assert q.CONTROL.ADSI() == task
    assert q._definition.executable == "kkrgen"
    assert not q._definition.mpi
    assert pickle.loads(pickle.dumps(q)).to_string() == q.to_string()
    if task == "PLOTSFN":
        assert q.SCF.FULLPOT()


def test_shape_functions_need_fullpot():
    p = InputParameters.create("PLOTSFN")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.SCF.FULLPOT.set(False)
    with pytest.raises((DataValidityError, DataValidityErrors), match="FULLPOT"):
        p.validate("save")
