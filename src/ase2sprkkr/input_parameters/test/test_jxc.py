import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("JXC")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    assert p.TASK.TASK() == "JXC"
    assert p.TASK.JXC() and not p.TASK.DMI()
    assert p.TAU.BZINT() == "POINTS" and p.ENERGY.NE().tolist() == [30]
    text = p.to_string(validate=True)
    assert "\n\tJXC\n" in text and "\n\tDMI\n" not in text
    assert "CLURAD=2.2" in text and "QREPONLY" not in text
    path = tmp_path / "jxc.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.JXC() and not q.TASK.DMI()
    assert q.to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_jxc_and_dmi_can_be_combined():
    p = parameters()
    p.TASK.DMI = True
    text = p.to_string(validate=True)
    assert "\n\tJXC\n" in text and "\n\tDMI\n" in text

    q = p._definition.read_from_string(text)
    assert q.TASK.JXC() and q.TASK.DMI()


def test_at_least_one_calculation_is_required():
    p = parameters()
    p.TASK.JXC = False
    with pytest.raises(DataValidityError, match="At least one of TASK.JXC or TASK.DMI"):
        p.validate("save")


def test_sp_srel_supported_for_jxc_only():
    p = parameters()
    p.MODE.MODE = "SP-SREL"
    p.validate("save")

    p.TASK.DMI = True
    with pytest.raises(DataValidityError, match="DMI requires magnetic FREL"):
        p.validate("save")


def test_cluster_options(tmp_path):
    p = parameters()
    p.TASK.set({"NSHLCLU": 2, "SELECT": ["Fe", "Co"], "QREPONLY": True})
    path = tmp_path / "selected.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.SELECT().tolist() == ["Fe", "Co"] and q.TASK.QREPONLY()
    assert q.TASK.NSHLCLU() == 2
    p.TASK.EXCLUDE = ["Ni"]
    with pytest.raises(DataValidityError, match="mutually exclusive"):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("GRID", 3), ("GRID", [5, 3]), ("NE", [32, 10]), ("SPLITSS", True)])
def test_unsupported_contour(name, value):
    p = parameters()
    p.ENERGY[name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("value", ["WEYL", "CLUSTER"])
def test_unsupported_integration(value):
    p = parameters()
    p.TAU.BZINT = value
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("value", [0, -1])
def test_invalid_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value


@pytest.mark.filterwarnings("ignore::ase2sprkkr.common.warnings.DataValidityWarning")
@pytest.mark.parametrize("section,name,value", [("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"), ("CONTROL", "NONMAG", True)])
def test_nonmagnetic_modes_rejected(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError, match="magnetic|magnetic calculation"):
        p.validate("save")
