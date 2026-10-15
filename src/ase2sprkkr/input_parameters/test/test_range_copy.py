import numpy as np
import pytest

from ...common.grammar_types import Range
from ..input_parameters import InputParameters


def test_range_comparison_handles_array_values():
    range_type = Range(float)

    assert range_type.is_the_same_value(
        np.array([45.0, 45.0]), np.array([45.0, 45.0])
    )
    assert not range_type.is_the_same_value(
        np.array([45.0, 45.0]), np.array([0.0, 0.0])
    )


@pytest.mark.parametrize("task", ["ARPES", "SPLEED"])
def test_copy_values_handles_explicit_range(task):
    parameters = InputParameters.create(task)
    parameters.SPEC_EL.THETA = [45.0, 45.0]

    copied = parameters.copy(copy_values=True)

    np.testing.assert_array_equal(copied.SPEC_EL.THETA(), [45.0, 45.0])
