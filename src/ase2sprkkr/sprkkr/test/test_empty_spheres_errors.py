from unittest.mock import Mock

from ase.build import bulk
import pytest

from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr.calculator import SPRKKR


def _calculator_without_potential():
    calculator = SPRKKR(atoms=bulk("Fe"), potential=False)
    calculator.input_parameters.CONTROL.POTFIL = "unused.pot"
    return calculator


@pytest.mark.parametrize(
    "empty_spheres",
    [True, {"min_radius": 0.5}, (0.5, 1.0)],
)
def test_explicit_empty_sphere_search_errors_propagate(monkeypatch, empty_spheres):
    calculator = _calculator_without_potential()

    def fail(*args, **kwargs):
        raise RuntimeError("empty-sphere search failed")

    monkeypatch.setattr(calculation_module, "add_empty_spheres", fail)

    with pytest.raises(RuntimeError, match="empty-sphere search failed"):
        calculator.save_input(potential=False, empty_spheres=empty_spheres)


def test_empty_spheres_false_skips_search(monkeypatch, tmp_path):
    calculator = _calculator_without_potential()
    search = Mock(side_effect=AssertionError("empty-sphere search must not run"))
    monkeypatch.setattr(calculation_module, "add_empty_spheres", search)

    calculator.save_input(
        potential=False,
        empty_spheres=False,
        directory=tmp_path,
    )

    search.assert_not_called()
