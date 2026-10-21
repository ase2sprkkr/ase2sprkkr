"""Regression for Python 3 arithmetic of structured result data."""

import numpy as np
import pytest
from pathlib import Path

from .. import ResultFile
from ..definitions.dos import DOS


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def bsf_result():
    return ResultFile.from_file(str(EXAMPLES / "SnTe_BLOCHSF.bsf"), unknown=False)


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


def test_left_scalar_multiplication_never_mutates_dos():
    original = _dos()
    before = original.dos.copy()
    scaled = 3 * original
    np.testing.assert_allclose(scaled.dos, 3 * before)
    np.testing.assert_array_equal(original.dos, before)
    assert not np.shares_memory(scaled.dos, original.dos)


def test_dos_rejects_different_energy_grids_for_all_binary_ops():
    import pytest

    left = _dos()
    right = DOS(left.energy + 0.25, left.dos.copy())
    for operation in (
        lambda: left + right,
        lambda: left - right,
        lambda: left * right,
        lambda: left / right,
        lambda: left.__iadd__(right),
        lambda: left.__imul__(right),
    ):
        with pytest.raises(ValueError, match="energy grids"):
            operation()
    np.testing.assert_array_equal(left.dos, [2.0, 4.0, 6.0])


def test_dos_rejects_mismatched_shapes():
    import pytest

    left = _dos()
    right = DOS(left.energy, np.ones((2, 3)))
    with pytest.raises(ValueError, match="different shapes"):
        left * right


def test_dos_scalar_add_and_subtract():
    dos = _dos()
    np.testing.assert_allclose((dos + 1).dos, [3.0, 5.0, 7.0])
    np.testing.assert_allclose((dos - 1).dos, [1.0, 3.0, 5.0])
    np.testing.assert_allclose(dos.dos, [2.0, 4.0, 6.0])


def test_bsf_arithmetic_preserves_source_and_checks_energy_grid(bsf_result):
    file = bsf_result
    original = file.RAW_DATA().copy()

    scaled = 2 * file
    np.testing.assert_allclose(scaled.RAW_DATA(), original * 2)
    np.testing.assert_array_equal(file.RAW_DATA(), original)

    repeated = file + file
    np.testing.assert_allclose(repeated.RAW_DATA(), original * 2)
    np.testing.assert_array_equal(file.RAW_DATA(), original)

    shifted = file.copy(copy_values=True)
    shifted["EMIN"].set(file.EMIN() + 0.01)
    with pytest.raises(ValueError, match="BSF E grids"):
        file * shifted
    np.testing.assert_array_equal(file.RAW_DATA(), original)


def test_bsf_rejects_different_spin_channels(bsf_result):
    normal = bsf_result
    spol = ResultFile.from_file(str(EXAMPLES / "SnTe_BLOCHSF_spol.bsf"), unknown=False)
    with pytest.raises(ValueError, match="BSF KEYWORD"):
        normal + spol


def test_arpes_arithmetic_preserves_geometry_and_recomputes_polarization():
    from pathlib import Path
    import pytest
    from .. import ResultFile

    examples = Path(__file__).resolve().parents[1] / "examples"
    file = ResultFile.from_file(
        str(examples / "Mo2Se4_ARPES_ARPES_data.spc"), unknown=False,
    )
    original = file.RAW_DATA().copy()

    for scaled in (2 * file, file * 2, file / 2, file + file):
        # The intensity channels change, but the original geometry does not.
        np.testing.assert_array_equal(scaled.ENERGY(), file.ENERGY())
        np.testing.assert_array_equal(scaled.THETA(), file.THETA())
        np.testing.assert_array_equal(scaled.K(), file.K())
        up, down = scaled.UP(), scaled.DOWN()
        expected = np.divide(
            up - down, up + down, out=np.zeros_like(up), where=(up + down) != 0,
        )
        np.testing.assert_allclose(scaled.POLARIZATION(), expected, atol=1e-12)
        np.testing.assert_array_equal(file.RAW_DATA(), original)

    shifted = file.copy(copy_values=True)
    shifted.RAW_DATA()[:, 1] += 0.1
    with pytest.raises(ValueError, match="ARPES ENERGY"):
        file * shifted
    np.testing.assert_array_equal(file.RAW_DATA(), original)


def test_dos_validation_has_single_override_for_scalar_array_and_grid():
    import pytest

    left = _dos()
    assert left._check_arithmetic(2) is None
    assert left._check_arithmetic(np.array([1.0, 1.0, 1.0])) is None
    assert left._check_arithmetic(_dos()) is None
    with pytest.raises(TypeError, match="Cannot combine"):
        left._check_arithmetic("invalid")
    with pytest.raises(ValueError, match="energy grids"):
        left._check_arithmetic(DOS(left.energy + 1, left.dos.copy()))


def test_result_arithmetic_does_not_use_legacy_assert_hook():
    from ..result_files import Arithmetic

    assert "_assert_arithmetic" not in Arithmetic.__dict__


def test_dos_rejects_array_broadcasting():
    import pytest

    dos = _dos()
    with pytest.raises(ValueError, match="different shapes"):
        dos *= np.ones((2, 3))
    np.testing.assert_array_equal(dos.dos, [2.0, 4.0, 6.0])


def test_bsf_grid_check_does_not_expect_nonexistent_shape_field(bsf_result):
    bsf = bsf_result
    assert bsf._check_arithmetic(bsf.copy(copy_values=True)) is None

def test_dos_numpy_array_arithmetic_uses_elementwise_values():
    original = _dos()
    factors = np.array([2.0, 3.0, 4.0])
    scaled = original * factors
    np.testing.assert_array_equal(scaled.dos, [4.0, 12.0, 24.0])
    np.testing.assert_array_equal(original.dos, [2.0, 4.0, 6.0])
    original *= factors
    np.testing.assert_array_equal(original.dos, [4.0, 12.0, 24.0])


def test_bsf_const_e_origin_is_part_of_arithmetic_grid():
    import pytest
    from ..definitions.bsf import BSFResultFile

    class _GridOnlyBSF(BSFResultFile):
        def __init__(self, start):
            self._grid = {
                "MODE": "CONST-E", "KEYWORD": "BSF", "NQ_EFF": 1,
                "NK1": 2, "NK2": 2, "VECK_START": start,
                "VECK1": np.array([1., 0., 0.]),
                "VECK2": np.array([0., 1., 0.]),
                "ERYD": np.array([0., 0.]),
            }

        def __getitem__(self, name):
            return lambda: self._grid[name]

        def RAW_DATA(self):
            return np.ones((4,))

        def MODE(self):
            return self._grid["MODE"]

    first = _GridOnlyBSF(np.array([0., 0., 0.]))
    shifted = _GridOnlyBSF(np.array([0.1, 0., 0.]))
    with pytest.raises(ValueError, match="BSF VECK_START grids"):
        first._check_arithmetic(shifted)
