import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("XAS")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_round_trip(tmp_path):
    p = parameters()
    assert not p.task_definition.mpi and p.task_definition.executable == "kkrgen"
    assert p.ENERGY.EMAX() == 4.0
    text = p.to_string(validate=True)
    assert "EMAX=4.0" in text and "FRAMETET=" not in text and "OUTPUT=" not in text
    assert "MECHECK" not in p.TASK
    path = tmp_path / "xas.in"
    p.save_to_file(path, atoms=False)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text


def test_explicit_options(tmp_path):
    p = parameters()
    p.TASK.set({"IT": 2, "CL": "3D", "POLOBS": "C", "OUTPUT": "SIGMA", "FRAMEGAM": 20.0})
    p.ENERGY.EMAX = 2.0
    path = tmp_path / "options.in"
    p.save_to_file(path, atoms=False)
    q = InputParameters.from_file(path)
    assert q.TASK.POLOBS() == "C" and q.TASK.OUTPUT() == "SIGMA"
    assert q.TASK.FRAMEGAM() == 20.0 and q.ENERGY.EMAX() == 2.0


@pytest.mark.parametrize("section,name,value", [
    ("MODE", "MODE", "NREL"), ("MODE", "MODE", "SREL"), ("MODE", "MODE", "SP-SREL"),
    ("ENERGY", "GRID", [6, 3]), ("ENERGY", "NE", [180, 30]), ("ENERGY", "SPLITSS", True),
    ("TASK", "CL", "1P"), ("TASK", "FRAMEGAM", float("inf")),
])
def test_invalid_combinations(section, name, value):
    p = parameters()
    p[section][name].set(value)
    with pytest.raises(DataValidityError):
        p.validate("save")


@pytest.mark.parametrize("name,value", [("IT", 0), ("CL", "6S"), ("DIPOL", 2)])
def test_invalid_values(name, value):
    with pytest.raises(DataValidityError):
        parameters().TASK[name].set(value)


def test_nab_override():
    p = parameters()
    p.TASK.ME = "NAB"
    with pytest.warns(DataValidityWarning, match="forces DIPOL"):
        p.validate("save")


def test_nondipole_precedence():
    p = parameters()
    p.TASK.NONDIPE2 = True
    p.TASK.NONDIPM1 = True
    with pytest.warns(DataValidityWarning, match="takes precedence"):
        p.validate("save")
