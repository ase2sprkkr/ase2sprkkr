import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("XES")
    p.CONTROL.POTFIL = "Cu.pot"
    return p


def test_round_trip(tmp_path):
    p = parameters()
    assert not p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    text = p.to_string(validate=True)
    assert "EMAX=" not in text and "POLOBS=" not in text and "DIPOL=" not in text
    path = tmp_path / "xes.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_options(tmp_path):
    p = parameters()
    p.TASK.set({"IT": 2, "CL": "3D", "POLOBS": "C", "TAUCORE": [.4, .6], "FRAMETET": 30., "WRREGIRR": True})
    path = tmp_path / "options.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.TAUCORE().tolist() == [.4, .6] and q.TASK.POLOBS() == "C"
    assert q.TASK.IT() == 2 and q.TASK.WRREGIRR()


@pytest.mark.parametrize("section,name,value", [
    ("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"), ("MODE", "MODE", "SP-SREL"),
    ("ENERGY", "GRID", [3, 5]), ("ENERGY", "NE", [300, 30]), ("ENERGY", "SPLITSS", True),
    ("TASK", "CL", "1P"), ("TASK", "TAUCORE", [-.1]),
])
def test_invalid_combinations(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("IT", 0), ("DIPOL", 2)])
def test_invalid_values(name, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set(value)


def test_nab_override():
    p = parameters()
    p.TASK.ME = "NAB"
    with pytest.warns(DataValidityWarning, match="forces DIPOL"):
        p.validate("save")


@pytest.mark.parametrize("value", [0, 1])
def test_energy_mesh(value):
    with pytest.raises(DataValidityError):
        parameters().ENERGY.NE = value


def test_core_shell_table_limit():
    with pytest.raises(DataValidityError):
        parameters().TASK.CL = "6S"
