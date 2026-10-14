import numpy as np
import pytest

from ..grammar_types import Array, AlwaysQString, Bool, Integer
from ..grammar_types.mixed import CustomMixed


def _type():
    return CustomMixed(
        (Bool.I, Array(AlwaysQString.I), Integer.I, AlwaysQString.I),
        string_type=AlwaysQString.I,
    )


def test_custom_mixed_dispatches_declared_python_value_types():
    value_type = _type()

    assert value_type.convert(True) is True
    assert value_type.convert(8) == 8
    assert value_type.convert("auto") == "auto"
    assert np.array_equal(value_type.convert(["srun", "-n", "8"]), ["srun", "-n", "8"])


def test_custom_mixed_rejects_unsupported_mappings():
    with pytest.raises(ValueError):
        _type().convert({"nonsense": True})
