"""Strict TASK-selector input autodetection, with a declarative cached IPD."""

import io
from pathlib import Path

import pyparsing as pp
import pytest

from ..input_parameters import InputParameters
from ..input_parameters_definitions import (
    InputParametersDefinition,
    MAX_INPUT_FILE_SIZE,
    _read_input_text,
    _task_preparse_definition,
    _task_selector_map,
)


_EXAMPLES = Path(__file__).parent.parent / "examples"


def test_preparse_ipd_is_cached_and_contains_only_repeated_task():
    a = _task_preparse_definition()
    assert a is _task_preparse_definition()
    assert a.grammar() is a.grammar()
    assert list(a.sections) == ["TASK"]
    assert list(a["TASK"].values) == ["TASK"]
    assert a["TASK"]["TASK"].is_repeated
    assert a["TASK"]["TASK"].type.keywords

    one = a.read_from_string("TASK ARPES")
    assert list(one.TASK.TASK()) == ["ARPES"]
    both = a.read_from_string("TASK DMI JXC")
    assert list(both.TASK.TASK()) == ["DMI", "JXC"]
    none = a.read_from_string("CONTROL ADSI=unrelated")
    assert not none.TASK.TASK()


def test_mapping_is_declarative_and_independent_of_adsi():
    mapping = _task_selector_map()
    assert mapping["JXC"] == mapping["DMI"]
    assert mapping["ARPES"] == mapping["SPLEED"]
    assert mapping["ARPES"] == mapping["AIPES"]
    assert mapping["ARPES"] == mapping["BAND"]
    assert mapping["BSF"] == (InputParameters.definition("BSF"),)
    for definition in InputParameters.definitions.values():
        assert definition.accepted_task_selectors()


def test_class_vs_instance_from_file():
    path = _EXAMPLES / "jxc.in"
    selected = InputParametersDefinition.from_file(path)
    explicit = InputParameters.definition("JXC").from_file(path)
    conventional = InputParameters.from_file(path)
    assert selected.definition_name == "JXC"
    assert conventional.definition_name == "JXC"
    assert explicit.definition_name == "JXC"
    assert selected.TASK.DMI()


def test_missing_task_selects_scf_and_not_adsi():
    path = _EXAMPLES / "scf.in"
    text = path.read_text()
    assert InputParameters.from_file(io.StringIO(text)).definition_name == "SCF"


def test_explicit_ipd_never_autodetects():
    path = _EXAMPLES / "jxc.in"
    scf = InputParameters.definition("SCF")
    with pytest.raises(pp.ParseBaseException):
        scf.from_file(path)


def test_limit_from_path_and_binary_and_text_stream(tmp_path):
    data = b"x" * (MAX_INPUT_FILE_SIZE + 1)
    path = tmp_path / "oversize.inp"
    path.write_bytes(data)
    for source in (path, io.BytesIO(data), io.StringIO("x" * (MAX_INPUT_FILE_SIZE + 1))):
        with pytest.raises(ValueError, match="exceeds"):
            InputParametersDefinition.from_file(source)
    assert _read_input_text(io.BytesIO(b"abc")) == "abc"


def test_unknown_selector_never_guesses():
    with pytest.raises(pp.ParseBaseException):
        InputParametersDefinition.from_file(io.StringIO("TASK UNKNOWN"))


def test_selected_parser_failure_never_falls_back():
    with pytest.raises(pp.ParseBaseException):
        InputParametersDefinition.from_file(io.StringIO("TASK ARPES INVALID"))


def test_error_policy_with_multiple_candidate_ipds(monkeypatch):
    from .. import input_parameters_definitions as mod

    class FakeDefinition:
        def __init__(self, name, outcome):
            self.name, self.outcome, self.calls = name, outcome, 0

        def read_from_string(self, text, allow_dangerous=False):
            self.calls += 1
            if isinstance(self.outcome, Exception):
                raise self.outcome
            return self.outcome

    first = FakeDefinition("FIRST", pp.ParseException("mismatch"))
    success = object()
    second = FakeDefinition("SECOND", success)
    unrelated = FakeDefinition("OTHER", object())
    monkeypatch.setattr(mod, "_task_selector_map", lambda: {
        "FIRST": (first,), "SECOND": (second,), "OTHER": (unrelated,),
    })

    assert mod._read_selected_definitions("content", ("FIRST", "SECOND")) is success
    assert (first.calls, second.calls, unrelated.calls) == (1, 1, 0)

    first.outcome = object()
    with pytest.raises(ValueError, match="Ambiguous"):
        mod._read_selected_definitions("content", ("FIRST", "SECOND"))

    first.outcome = TypeError("programmer error")
    with pytest.raises(TypeError, match="programmer error"):
        mod._read_selected_definitions("content", ("FIRST", "SECOND"))
