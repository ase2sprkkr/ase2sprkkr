"""Regression for Python 3 arithmetic of structured result data."""

import numpy as np

from ..definitions.dos import DOS


def _dos():
    return DOS(
        energy=np.array([-1.0, 0.0, 1.0]),
        dos=np.array([2.0, 4.0, 6.0]),
    )


def test_dos_true_division_returns_independent_values():
    original = _dos()
    divided = original / 2.0
    np.testing.assert_array_equal(divided.dos, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(original.dos, [2.0, 4.0, 6.0])


def test_dos_in_place_true_division_updates_values():
    original = _dos()
    values = original.dos
    original /= 2.0
    assert original.dos is values
    np.testing.assert_array_equal(original.dos, [1.0, 2.0, 3.0])


def test_dos_true_division_by_dos_uses_numerator_values():
    left = _dos()
    right = DOS(left.energy, np.array([2.0, 2.0, 3.0]))
    divided = left / right
    np.testing.assert_array_equal(divided.dos, [1.0, 2.0, 2.0])
    np.testing.assert_array_equal(left.dos, [2.0, 4.0, 6.0])
