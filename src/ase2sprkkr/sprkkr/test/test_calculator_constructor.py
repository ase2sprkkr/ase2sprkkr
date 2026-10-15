import pytest

from ...input_parameters.input_parameters import InputParameters
from ...potentials.potentials import Potential
from ..calculator import SPRKKR


def test_constructor_kwargs_configure_input_parameters():
    calculator = SPRKKR(NITER=7)

    assert calculator.get("NITER") == 7


def test_constructor_kwargs_override_options():
    calculator = SPRKKR(options={"NITER": 3}, NITER=5)

    assert calculator.get("NITER") == 5


def test_constructor_rejects_unknown_parameter_keyword():
    with pytest.raises(KeyError, match="NOT_A_SPRKKR_OPTION"):
        SPRKKR(NOT_A_SPRKKR_OPTION=1)


def test_calculator_convenience_aliases_are_always_available():
    assert SPRKKR.InputParameters is InputParameters
    assert SPRKKR.Potential is Potential
