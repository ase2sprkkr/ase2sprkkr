"""PLOT2DPR window syntax and native division-by-zero hazards."""
import pickle
import pytest
from ase2sprkkr.common.warnings import DataValidityError, DataValidityWarning
from ase2sprkkr.input_parameters.input_parameters import InputParameters


def parameters():
    p = InputParameters.create("PLOT2DPR")
    p.CONTROL.POTFIL = "Fe.pot"
    return p


def test_defaults_and_file_round_trip(tmp_path):
    p = parameters()
    assert p._definition.executable == "kkrgen" and not p.is_mpi()
    assert p.TASK.AVEC().tolist() == [2.0, 0.0, 0.0]
    assert p.TASK.BVEC().tolist() == [0, 2, 0]
    assert p.TASK.DVEC().tolist() == [0.0, 0.0, 0.0]
    assert p.TASK.NA() == 50 and p.TASK.NB() is None
    assert p.ENERGY.NE().tolist() == [1]
    for name in ("AVEC", "BVEC", "DVEC", "NA", "NB"):
        assert f"\t{name}=" not in p.TASK.to_string()
    path = tmp_path / "plot2dpr.in"
    p.save_to_file(path, atoms=False)
    text = p.to_string(validate=True)
    assert InputParameters.from_file(path).to_string() == text
    assert pickle.loads(pickle.dumps(p)).to_string() == text
    assert p.copy(copy_values=True).to_string() == text


def test_custom_window_round_trip():
    p = parameters()
    p.TASK.set({"AVEC": [3.0, 0.0, 0.0], "BVEC": [0.0, 1.0, 0.0], "DVEC": [-1.0, 0.0, 0.0], "NA": 60, "NB": 20})
    text = p.to_string(validate=True)
    assert "NB=20" in text
    assert p._definition.read_from_string(text).to_string() == text


@pytest.mark.parametrize("keyword", ["PLOT2DPR", "PLOT2DPOT"])
def test_plot_selection(keyword, tmp_path):
    p = parameters()
    p.TASK.TASK = keyword
    text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).TASK.TASK() == keyword
    path = tmp_path / "selection.in"
    path.write_text(text)
    assert InputParameters.from_file(path).TASK.TASK() == keyword


def test_density_selection_warns_about_asa():
    p = parameters()
    p.TASK.TASK = "PLOT2DRHO"
    with pytest.warns(DataValidityWarning, match="still writes potential"):
        text = p.to_string(validate=True)
    assert p._definition.read_from_string(text).TASK.TASK() == "PLOT2DRHO"


@pytest.mark.parametrize("a,b", [([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), ([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
                                 ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])])
def test_degenerate_window_rejected(a, b):
    p = parameters()
    p.TASK.set({"AVEC": a, "BVEC": b})
    with pytest.raises(DataValidityError, match="nonzero, nonparallel"):
        p.validate("save")


def test_orthogonalization_warns():
    p = parameters()
    p.TASK.BVEC = [1.0, 2.0, 0.0]
    with pytest.warns(DataValidityWarning, match="orthogonalizes"):
        p.validate("save")


def test_auto_nb_must_not_round_to_zero():
    p = parameters()
    p.TASK.BVEC = [0.0, 0.001, 0.0]
    with pytest.raises(DataValidityError, match="rounds to zero"):
        p.validate("save")
    p.TASK.NB = 1
    p.validate("save")


@pytest.mark.parametrize("key", ["NA", "NB"])
def test_zero_intervals_rejected(key):
    with pytest.raises(DataValidityError):
        parameters().TASK[key].set(0)


@pytest.mark.parametrize("op", ["DMFT", "LDA+U"])
def test_correlations_rejected(op):
    p = parameters()
    p.MODE.OP = op
    with pytest.raises(DataValidityError, match="does not support"):
        p.validate("save")


@pytest.mark.parametrize("key", ["AVEC", "BVEC"])
def test_window_vectors_require_three_components(key):
    with pytest.raises(DataValidityError):
        parameters().TASK[key].set([1.0, 2.0])
