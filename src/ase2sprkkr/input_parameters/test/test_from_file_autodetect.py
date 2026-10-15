import pyparsing as pp
import pytest

from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError


class _Definition:
    def __init__(self, outcome):
        self.outcome = outcome
        self.calls = 0

    def read_from_file(self, filename, allow_dangerous=False):
        self.calls += 1
        if isinstance(self.outcome, BaseException):
            raise self.outcome
        return self.outcome


def _parameters_with(definitions):
    class Parameters(InputParameters):
        pass

    Parameters.definitions = definitions
    return Parameters


def test_from_file_tries_next_definition_after_parse_mismatch():
    mismatch = _Definition(pp.ParseException("not this task"))
    expected = object()
    match = _Definition(expected)
    parameters = _parameters_with({"first": mismatch, "second": match})

    assert parameters.from_file("input.inp") is expected
    assert mismatch.calls == 1
    assert match.calls == 1


def test_from_file_tries_next_definition_after_data_validity_error():
    mismatch = _Definition(DataValidityError("invalid parsed data"))
    expected = object()
    match = _Definition(expected)
    parameters = _parameters_with({"first": mismatch, "second": match})

    assert parameters.from_file("input.inp") is expected
    assert mismatch.calls == 1
    assert match.calls == 1


@pytest.mark.parametrize("error", [
    AttributeError("parser bug"),
    TypeError("parser bug"),
    ValueError("parser bug"),
    KeyError("parser bug"),
])
def test_from_file_does_not_mask_programming_errors(error):
    broken = _Definition(error)
    fallback = _Definition(object())
    parameters = _parameters_with({"broken": broken, "fallback": fallback})

    with pytest.raises(type(error), match="parser bug"):
        parameters.from_file("input.inp")

    assert broken.calls == 1
    assert fallback.calls == 0
