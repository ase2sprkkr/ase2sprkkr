from unittest.mock import Mock

from ase.build import bulk
import pytest

from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr.calculator import SPRKKR


def _calculator_without_potential():
    calculator = SPRKKR(atoms=bulk("Fe"), potential=False)
    calculator.input_parameters.CONTROL.POTFIL = "unused.pot"
    return calculator


def test_explicit_empty_sphere_search_errors_propagate(monkeypatch):
    # All failures occur before any input files are written. Reusing the
    # calculator also verifies preparation remains usable after each error.
    calculator = _calculator_without_potential()
    search = Mock(side_effect=RuntimeError("empty-sphere search failed"))
    monkeypatch.setattr(calculation_module, "add_empty_spheres", search)

    for index, (empty_spheres, expected_options) in enumerate(
        (
            (True, {}),
            ({"min_radius": 0.5}, {"min_radius": 0.5}),
            ((0.5, 1.0), {"min_radius": 0.5, "max_radius": 1.0}),
        ),
        start=1,
    ):
        with pytest.raises(RuntimeError, match="empty-sphere search failed"):
            calculator.save_input(potential=False, empty_spheres=empty_spheres)

        assert search.call_count == index
        assert search.call_args.kwargs == expected_options


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
