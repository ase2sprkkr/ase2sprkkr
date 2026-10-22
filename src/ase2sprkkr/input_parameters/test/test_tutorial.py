import io
import pytest
from ..input_parameters import InputParameters
from ..definitions.tutorial import TUTORIAL_TASKS
from ...common.warnings import DataValidityError, DataValidityErrors


def parameters(task="YLM"):
    p = InputParameters.create("TUTORIAL")
    p.CONTROL.POTFIL.set("Fe.pot")
    p.TASK["TUT-TASK"].set(task)
    if task == "DOS-CLU":
        p.TAU.CLUSTER.set(True)
        p.TAU.CLURAD.set(2.2)
    return p


@pytest.mark.parametrize("task", list(TUTORIAL_TASKS))
def test_roundtrip(task):
    p = parameters(task)
    q = p._definition.read_from_file(io.StringIO(p.to_string(validate=True)))
    assert q.TASK["TUT-TASK"]() == task
    assert q.task_definition.executable == "kkrgen"
    assert not q.task_definition.mpi


@pytest.mark.parametrize("task,option,value", [("YLM", "M", 3), ("YLM", "L", -1), ("JLX", "NX", 1), ("JLX", "LMAX", 5), ("JLX", "XMAX", .001), ("JLX", "XMIN", float("nan")), ("JLX-ASYMPT", "XMIN", .1)])
def test_invalid(task, option, value):
    p = parameters(task)
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.TASK[option].set(value)
        p.validate("save")


@pytest.mark.parametrize("task", ["POT-DECOMP", "POT-XC", "CORE-ALL"])
def test_potential_decomposition_rejects_fullpot(task):
    p = parameters(task)
    p.SCF.FULLPOT.set(True)
    with pytest.raises((DataValidityError, DataValidityErrors), match="ASA"):
        p.validate("save")


@pytest.mark.parametrize("task", ["POT-DECOMP", "POT-XC", "CORE-ALL"])
def test_potential_decomposition_ignores_mathematical_range(task):
    p = parameters(task)
    p.TASK.M.set(3)
    p.TASK.XMAX.set(.001)
    p.validate("save")


@pytest.mark.parametrize("task", ["RWF-L", "RWF-E"])
def test_rwf_defaults_and_pickle(task):
    import pickle
    p = parameters(task)
    assert p.MODE.MODE() == "SREL"
    assert p.ENERGY.NE().tolist() == [10]
    assert p.ENERGY.EMIN() == 1e-7
    q = pickle.loads(pickle.dumps(p))
    assert q.to_string(validate=True) == p.to_string(validate=True)


@pytest.mark.parametrize("task,option,value", [
    ("RWF-E", "ENERGY.NE", [1]), ("RWF-L", "ENERGY.EMINEV", 1.),
    ("RWF-E", "MODE.MODE", "FREL"), ("RWF-L", "SCF.FULLPOT", True),
    ("RWF-E", "TASK.ITSEL", 0), ("DOS-BULK", "TAU.CLUSTER", True),
    ("DOS-SS", "TAU.KKRMODE", "TB"), ("DOS-FEG", "TAU.CLUSTER", True),
    ("DOS-BULK", "ENERGY.NE", [1]),
])
def test_remaining_branch_validation(task, option, value):
    p = parameters(task)
    with pytest.raises((DataValidityError, DataValidityErrors, ValueError)):
        p.set({option: value})
        p.validate("save")


def test_cluster_name_does_not_configure_cluster():
    p = parameters("DOS-BULK")
    p.TASK["TUT-TASK"].set("DOS-CLU")
    with pytest.raises((DataValidityError, DataValidityErrors), match="DOS-CLU"):
        p.validate("save")


def test_changing_subtask_updates_written_defaults():
    p = parameters()
    assert p.ENERGY.NE().tolist() == [1]
    p.TASK["TUT-TASK"].set("DOS-SS")
    assert p.ENERGY.NE().tolist() == [300]
    assert p.TAU.KKRMODE() == "STANDARD"
    p.TASK["TUT-TASK"].set("RWF-E")
    assert p.ENERGY.NE().tolist() == [10]
    assert p.MODE.MODE() == "SREL"
