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
    assert second.calls == 0


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


def _arpes_input_text():
    parameters = InputParameters.create("ARPES")
    parameters.CONTROL.POTFIL.set("Fe.pot")
    return parameters.to_string(validate=False)


def test_task_section_rejects_another_registered_task_as_custom_flag():
    text = _arpes_input_text().replace(
        "\nTASK\n\tARPES\n",
        "\nTASK\n\tBSF\n",
        1,
    )

    with pytest.raises(pp.ParseBaseException):
        InputParameters.definition("ARPES").read_from_string(text)


def test_task_section_still_allows_user_custom_flags():
    text = _arpes_input_text().replace(
        "\nTASK\n\tARPES\n",
        "\nTASK\n\tARPES\n\tMY_CUSTOM_FLAG\n",
        1,
    )

    parsed = InputParameters.definition("ARPES").read_from_string(text)

    assert parsed.TASK["MY_CUSTOM_FLAG"]() is True


def test_from_file_detects_bsf_from_task_selector(tmp_path):
    parameters = InputParameters.create("BSFEK")
    parameters.CONTROL.POTFIL.set("Fe.pot")
    parameters.TASK.KPATH.set("6")
    filename = tmp_path / "bsf.inp"
    filename.write_text(parameters.to_string(validate=False))

    parsed = InputParameters.from_file(filename)

    assert parsed.task_name == "BSFEK"
    assert parsed.TASK.KPATH() == "6"
