from types import SimpleNamespace

import pyparsing as pp
import pytest

from ..input_parameters import InputParameters
from ...common.warnings import DataValidityError


class _Definition:
    def __init__(self, outcome, *, schema_adsi=None):
        self.outcome = outcome
        self.schema_adsi = schema_adsi
        self.calls = 0

    def __getitem__(self, name):
        assert name == "CONTROL"
        return {"ADSI": SimpleNamespace(default_value=self.schema_adsi)}

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


def _candidate(adsi, task_name):
    return SimpleNamespace(
        CONTROL=SimpleNamespace(ADSI=lambda: adsi),
        task_name=task_name,
    )


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


def test_from_file_prefers_success_matching_registered_adsi_task():
    permissive = _Definition(_candidate("SCF", "ARPES"))
    expected = _candidate("SCF", "SCF")
    exact = _Definition(expected)
    parameters = _parameters_with({"first": permissive, "second": exact})

    assert parameters.from_file("input.inp") is expected
    assert permissive.calls == 1
    assert exact.calls == 1


def test_from_file_keeps_first_success_when_adsi_is_not_a_task():
    expected = _candidate("custom label", "ARPES")
    first = _Definition(expected)
    second = _Definition(_candidate("another label", "SCF"))
    parameters = _parameters_with({"first": first, "second": second})

    assert parameters.from_file("input.inp") is expected
    assert first.calls == 1
    # Unknown ADSI can still be a schema discriminator; check all schemas
    # before falling back to the first successfully parsed input.
    assert second.calls == 1


def test_from_file_prefers_schema_adsi_for_shared_task_family():
    first = _Definition(_candidate("BLOCHSF", "ARPES"), schema_adsi="ARPES")
    expected = _candidate("BLOCHSF", "BSFEK")
    matching = _Definition(expected, schema_adsi="BLOCHSF")
    parameters = _parameters_with({"arpes": first, "bsf": matching})

    assert parameters.from_file("input.inp") is expected
    assert first.calls == 1
    assert matching.calls == 1


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
