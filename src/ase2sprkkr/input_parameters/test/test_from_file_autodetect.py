"""Input parsing selects one full IPD from TASK.TASK, defaulting to SCF."""

import io

import pyparsing as pp
import pytest

from ..input_parameters import InputParameters, _task_selector_definition


class _Definition:
    def __init__(self, selectors, result):
        self.selectors = selectors
        self.name = selectors[0]
        self.result = result
        self.calls = 0

    def task_selector_names(self):
        return self.selectors

    def read_from_file(self, filename, allow_dangerous=False):
        self.calls += 1
        if isinstance(self.result, Exception):
            raise self.result
        return self.result


def _parameters_with(definitions):
    class Parameters(InputParameters):
        pass

    Parameters.definitions = definitions
    return Parameters


def test_minimal_selector_ipd_is_cached_and_contains_only_task_task():
    a = _task_selector_definition(InputParameters)
    assert a is _task_selector_definition(InputParameters)
    assert list(a.sections) == ["TASK"]
    assert list(a["TASK"].values) == ["TASK"]

    parsed = a.read_from_string("TASK ARPES")
    assert tuple(parsed.TASK.TASK()) == ("ARPES",)
    assert tuple(a.read_from_string("TASK DMI JXC").TASK.TASK()) == ("DMI", "JXC")
    assert tuple(a.read_from_string("TASK JXC DMI").TASK.TASK()) == ("JXC", "DMI")
    assert a.read_from_string("CONTROL ADSI=SCF").TASK.TASK() is None


def test_from_file_uses_only_selected_input_definition():
    expected = object()
    wrong = _Definition(["DOS"], object())
    correct = _Definition(["ARPES", "SPLEED"], expected)
    parameters = _parameters_with({"wrong": wrong, "right": correct})
    parameters.task_registry = object()

    assert parameters.from_file(io.StringIO("TASK ARPES")) is expected
    assert wrong.calls == 0
    assert correct.calls == 1


def test_from_file_with_preceding_task_flag_selects_jxc():
    expected = object()
    jxc = _Definition(["DMI", "JXC"], expected)
    scf = _Definition(["SCF"], object())
    parameters = _parameters_with({"scf": scf, "jxc": jxc})

    assert parameters.from_file(io.StringIO("TASK DMI JXC")) is expected
    assert jxc.calls == 1
    assert scf.calls == 0


def test_from_file_without_task_selects_scf_only():
    expected = object()
    scf = _Definition(["SCF"], expected)
    other = _Definition(["DOS"], object())
    parameters = _parameters_with({"scf": scf, "other": other})

    assert parameters.from_file(io.StringIO("CONTROL ADSI=DOS")) is expected
    assert scf.calls == 1
    assert other.calls == 0


def test_from_file_propagates_selected_parser_error_without_fallback():
    failure = pp.ParseException("invalid selected schema")
    chosen = _Definition(["ARPES"], failure)
    other = _Definition(["DOS"], object())
    parameters = _parameters_with({"chosen": chosen, "other": other})

    with pytest.raises(pp.ParseException, match="invalid selected schema"):
        parameters.from_file(io.StringIO("TASK ARPES"))
    assert chosen.calls == 1
    assert other.calls == 0


def test_from_file_unknown_task_never_guesses():
    candidate = _Definition(["ARPES"], object())
    parameters = _parameters_with({"candidate": candidate})

    with pytest.raises(ValueError, match="input-parameter definition"):
        parameters.from_file(io.StringIO("TASK UNKNOWN"))
    assert candidate.calls == 0


def test_from_file_rejects_ambiguous_selector():
    first = _Definition(["ARPES"], object())
    second = _Definition(["ARPES"], object())
    parameters = _parameters_with({"first": first, "second": second})

    with pytest.raises(ValueError, match="multiple input parameter definitions"):
        parameters.from_file(io.StringIO("TASK ARPES"))
    assert first.calls == second.calls == 0


def test_from_file_ignores_control_adsi():
    result = object()
    selected = _Definition(["ARPES"], result)
    other = _Definition(["DOS"], object())
    parameters = _parameters_with({"selected": selected, "other": other})

    assert parameters.from_file(io.StringIO(
        "CONTROL ADSI=DOS\nTASK ARPES"
    )) is result
    assert selected.calls == 1
    assert other.calls == 0


def test_multiple_distinct_schema_candidates_try_only_the_selected_ipds():
    first = _Definition(["DMI"], pp.ParseException("wrong definition"))
    result = object()
    second = _Definition(["JXC"], result)
    unused = _Definition(["DOS"], object())
    cls = _parameters_with({"a": first, "b": second, "unused": unused})

    assert cls.from_file(io.StringIO("TASK DMI JXC")) is result
    assert (first.calls, second.calls, unused.calls) == (1, 1, 0)


def test_multiple_distinct_schema_candidates_reject_ambiguity():
    first = _Definition(["DMI"], object())
    second = _Definition(["JXC"], object())
    cls = _parameters_with({"a": first, "b": second})
    with pytest.raises(ValueError, match="Ambiguous TASK"):
        cls.from_file(io.StringIO("TASK DMI JXC"))
    assert (first.calls, second.calls) == (1, 1)


def test_multiple_distinct_schema_candidates_do_not_hide_programming_error():
    first = _Definition(["DMI"], TypeError("programming bug"))
    second = _Definition(["JXC"], object())
    cls = _parameters_with({"a": first, "b": second})
    with pytest.raises(TypeError, match="programming bug"):
        cls.from_file(io.StringIO("TASK DMI JXC"))
    assert (first.calls, second.calls) == (1, 0)


def test_multiple_distinct_schema_candidates_no_other_fallback():
    first = _Definition(["DMI"], pp.ParseException("one"))
    second = _Definition(["JXC"], pp.ParseException("two"))
    unused = _Definition(["DOS"], object())
    cls = _parameters_with({"a": first, "b": second, "unused": unused})
    with pytest.raises(ValueError, match="No selected IPD could parse"):
        cls.from_file(io.StringIO("TASK DMI JXC"))
    assert (first.calls, second.calls, unused.calls) == (1, 1, 0)
