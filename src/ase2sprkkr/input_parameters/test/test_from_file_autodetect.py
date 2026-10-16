from pathlib import Path
from types import SimpleNamespace

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


class _ExplicitParsed:
    task_definition = SimpleNamespace(
        name="TEST",
        input_parameters_definition_name="TEST",
    )

    def __getitem__(self, name):
        if name != "TASK":
            raise KeyError(name)
        return self

    def is_set(self):
        return True

    def __call__(self):
        return "TEST"


def _parameters_with(definitions):
    class Parameters(InputParameters):
        pass

    Parameters.definitions = definitions
    return Parameters


def test_from_file_tries_next_definition_after_parse_mismatch():
    mismatch = _Definition(pp.ParseException("not this task"))
    expected = _ExplicitParsed()
    match = _Definition(expected)
    parameters = _parameters_with({"first": mismatch, "second": match})

    assert parameters.from_file("input.inp") is expected
    assert mismatch.calls == 1
    assert match.calls == 1


def test_from_file_tries_next_definition_after_data_validity_error():
    mismatch = _Definition(DataValidityError("invalid parsed data"))
    expected = _ExplicitParsed()
    match = _Definition(expected)
    parameters = _parameters_with({"first": mismatch, "second": match})

    assert parameters.from_file("input.inp") is expected
    assert mismatch.calls == 1
    assert match.calls == 1


def test_from_file_requires_explicit_task_discriminator():
    filename = Path(__file__).parent.parent / "examples" / "jxc.in"

    parameters = InputParameters.from_file(filename)

    assert parameters.definition_name == "JXC"
    assert parameters.task_name == "JXC"
    assert parameters.TASK.TASK.is_set()


def test_read_for_task_rejects_file_with_foreign_task_identity():
    filename = Path(__file__).parent.parent / "examples" / "jxc.in"
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["DOS"]

    with pytest.raises(ValueError, match="Input file read for task 'DOS'"):
        registry.read_for_task(task, filename)


@pytest.mark.parametrize("error", [
    AttributeError("parser bug"),
    TypeError("parser bug"),
    ValueError("parser bug"),
    KeyError("parser bug"),
])
def test_from_file_does_not_mask_programming_errors(error):
    broken = _Definition(error)
    fallback = _Definition(_ExplicitParsed())
    parameters = _parameters_with({"broken": broken, "fallback": fallback})

    with pytest.raises(type(error), match="parser bug"):
        parameters.from_file("input.inp")

    assert broken.calls == 1
    assert fallback.calls == 0
