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
    _TaskSelectorMap,
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
    assert mapping.definitions("JXC", user_defined=False) == mapping.definitions(
        "DMI", user_defined=False
    )
    for task in ("ARPES", "SPLEED", "AIPES", "BAND", "BSF"):
        assert mapping.definitions(task, user_defined=False) == (
            InputParameters.definition(task),
        )
        assert mapping.definitions(task) == mapping[task]
        assert all(
            definition.is_user_defined
            for definition in mapping.definitions(task, user_defined=True)
        )
    for definition in InputParameters.definitions.values():
        assert definition.accepted_task_selectors()


def test_ipd_origin_defaults_registration_and_copy(tmp_path):
    from types import SimpleNamespace

    builtin = InputParameters.definition("ARPES")
    assert builtin.is_user_defined is False
    assert builtin.copy().is_user_defined is False

    custom = InputParametersDefinition("MYARPES", [])
    assert custom.is_user_defined is True
    assert custom.copy().is_user_defined is True

    module = SimpleNamespace(
        __name__="_user_myarpes",
        __file__=str(tmp_path / "myarpes.py"),
        __doc__="",
    )
    InputParameters.input_parameters_registry.prepare_value(
        "MYARPES", module, custom
    )
    assert custom.is_user_defined is True
    assert custom.copy().is_user_defined is True

    mapping = _TaskSelectorMap({"ARPES": (builtin, custom)})
    assert mapping["ARPES"] == (builtin, custom)
    assert mapping.definitions("ARPES") == (builtin, custom)
    assert mapping.definitions("ARPES", user_defined=False) == (builtin,)
    assert mapping.definitions("ARPES", user_defined=True) == (custom,)
    with pytest.raises(TypeError, match="user_defined"):
        mapping.definitions("ARPES", user_defined="False")


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


def test_explicit_ipd_never_autodetects(monkeypatch):
    from unittest.mock import Mock

    scf = InputParameters.definition("SCF")
    parse = Mock(return_value=object())
    monkeypatch.setattr(scf, "read_from_file", parse)
    expected = scf.from_file(io.StringIO("TASK ARPES"))
    assert expected is parse.return_value
    parse.assert_called_once()
    assert parse.call_args.args[0].read() == "TASK ARPES"


@pytest.fixture
def small_input_limit(monkeypatch):
    # Preserve the configured production cap but use a tiny boundary for
    # byte-counting cases that never need a megabyte of payload.
    from .. import input_parameters_definitions as definitions

    assert definitions.MAX_INPUT_FILE_SIZE == MAX_INPUT_FILE_SIZE
    assert MAX_INPUT_FILE_SIZE == 1024 * 1024
    limit = 64
    monkeypatch.setattr(definitions, "MAX_INPUT_FILE_SIZE", limit)
    return limit


def test_limit_from_path_and_binary_and_text_stream(tmp_path, small_input_limit):
    data = b"x" * (small_input_limit + 1)
    path = tmp_path / "oversize.inp"
    path.write_bytes(data)
    for source in (path, io.BytesIO(data), io.StringIO("x" * (small_input_limit + 1))):
        with pytest.raises(ValueError, match="exceeds"):
            InputParametersDefinition.from_file(source)
    assert _read_input_text(io.BytesIO(b"abc")) == "abc"


def test_unrecognized_task_switch_without_known_selector_defaults_to_scf():
    # Native TASK sections may also contain non-selector flags. When no
    # registered TASK selector occurs, the selected IPD must still be SCF.
    preview = _task_preparse_definition().read_from_string("TASK UNKNOWN")
    value = preview.TASK.TASK()
    assert value is None or len(value) == 0


def test_selected_parser_failure_never_falls_back(monkeypatch):
    from unittest.mock import Mock
    from .. import input_parameters_definitions as mod

    selected = InputParameters.definition("ARPES")
    fallback = InputParameters.definition("SCF")
    error = pp.ParseException("invalid ARPES input")
    chosen_reader = Mock(side_effect=error)
    unrelated_reader = Mock(return_value=object())

    monkeypatch.setattr(selected, "read_from_file", chosen_reader)
    monkeypatch.setattr(fallback, "read_from_file", unrelated_reader)
    monkeypatch.setattr(mod, "_task_selector_map", lambda: {
        "ARPES": (selected,),
        "SCF": (fallback,),
    })

    with pytest.raises(pp.ParseException, match="invalid ARPES"):
        InputParametersDefinition.from_file(io.StringIO("TASK ARPES"))

    chosen_reader.assert_called_once()
    unrelated_reader.assert_not_called()

def test_error_policy_with_multiple_candidate_ipds(monkeypatch):
    from .. import input_parameters_definitions as mod

    class FakeDefinition:
        def __init__(self, name, outcome):
            self.name, self.outcome, self.calls = name, outcome, 0

        def read_from_file(self, stream, allow_dangerous=False, *, filename=None):
            self.calls += 1
            stream.read()
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


def test_unordered_preparse_recognizes_interleaved_selectors():
    preview = _task_preparse_definition()
    cases = [
        "TASK DMI CLURAD=4.0 JXC",
        "TASK JXC CLURAD=4.0 DMI",
        "TASK DMI CLURAD=4.0 JXC NSHLCLU=3 DMI",
    ]
    for text in cases:
        parsed = preview.read_from_string(text)
        selected = list(parsed.TASK.TASK())
        assert selected == [
            token for token in text.split()
            if token in ("DMI", "JXC")
        ]


def test_class_autodetection_tries_all_interleaved_candidates(monkeypatch):
    from .. import input_parameters_definitions as mod

    class Candidate:
        def __init__(self, name, result):
            self.name = name
            self.result = result
            self.calls = []

        def read_from_file(self, stream, allow_dangerous=False, *, filename=None):
            self.calls.append((stream.read(), allow_dangerous))
            if isinstance(self.result, Exception):
                raise self.result
            return self.result

    first = Candidate("DMI-IPD", pp.ParseException("not a DMI file"))
    expected = object()
    second = Candidate("JXC-IPD", expected)
    unrelated = Candidate("SCF-IPD", object())
    monkeypatch.setattr(mod, "_task_selector_map", lambda: {
        "DMI": (first,), "JXC": (second,), "SCF": (unrelated,),
    })

    value = InputParametersDefinition.from_file(
        io.StringIO("TASK DMI CLURAD=4.0 JXC")
    )
    assert value is expected
    assert len(first.calls) == len(second.calls) == 1
    assert unrelated.calls == []
    assert first.calls[0][0] == second.calls[0][0]


def test_all_candidate_failures_keep_each_diagnostic(monkeypatch):
    from .. import input_parameters_definitions as mod

    class Failed:
        def __init__(self, name, error):
            self.name, self.error = name, error

        def read_from_file(self, stream, allow_dangerous=False, *, filename=None):
            raise pp.ParseException(self.error)

    monkeypatch.setattr(mod, "_task_selector_map", lambda: {
        "DMI": (Failed("IPD-DMI", "missing parameter DMI"),),
        "JXC": (Failed("IPD-JXC", "bad ENERGY grid"),),
    })
    with pytest.raises(ValueError) as exc:
        InputParametersDefinition.from_file(
            io.StringIO("TASK DMI CLURAD=4.0 JXC")
        )
    message = str(exc.value)
    assert "IPD-DMI: ParseException: missing parameter DMI" in message
    assert "IPD-JXC: ParseException: bad ENERGY grid" in message
    assert exc.value.__cause__ is not None


def test_named_source_preserved_for_class_and_instance():
    source = _EXAMPLES / "scf.in"
    for parser in (InputParametersDefinition, InputParameters.definition("SCF")):
        parsed = parser.from_file(source)
        assert parsed.original_filename == str(source)

    named = io.StringIO(source.read_text())
    named.name = "origin-from-stream.inp"
    parsed = InputParametersDefinition.from_file(named)
    assert parsed.original_filename == named.name


def test_binary_and_text_streams_are_read_only_once():
    class CountedStream(io.StringIO):
        def __init__(self, data):
            super().__init__(data)
            self.calls = 0

        def read(self, *args):
            self.calls += 1
            return super().read(*args)

    source = (_EXAMPLES / "scf.in").read_text()
    stream = CountedStream(source)
    parsed = InputParametersDefinition.from_file(stream)
    assert parsed.definition_name == "SCF"
    assert stream.calls == 1


def test_size_limit_exact_boundary_and_multibyte_utf8(small_input_limit):
    limit = small_input_limit
    assert len(_read_input_text(io.BytesIO(b"x" * limit))) == limit
    with pytest.raises(ValueError, match="exceeds"):
        _read_input_text(io.StringIO("č" * (limit // 2 + 1)))
    with pytest.raises(ValueError, match="exceeds"):
        _read_input_text(io.BytesIO("č".encode("utf-8") * (limit // 2 + 1)))


def test_preparse_uses_strict_grammar_when_full_parse_is_dangerous(monkeypatch):
    from .. import input_parameters_definitions as mod

    preparse = _task_preparse_definition()
    calls = []
    original_read = preparse.read_from_string

    def read(*args, **kwargs):
        calls.append(kwargs["allow_dangerous"])
        return original_read(*args, **kwargs)

    monkeypatch.setattr(preparse, "read_from_string", read)
    source = _EXAMPLES / "scf.in"
    InputParametersDefinition.from_file(source, allow_dangerous=True)
    assert calls == [False]


def test_scalar_task_selector_and_registration_immutability_contract():
    from .. import input_parameters_definitions as definitions

    schema = InputParametersDefinition("special", [])
    assert schema.accepted_task_selectors() == ("SPECIAL",)
    schema.task_selectors = "dmi"
    assert schema.accepted_task_selectors() == ("DMI",)
    schema.task_selectors = ("jxc", "dmi")
    assert schema.accepted_task_selectors() == ("JXC", "DMI")
    for invalid in ((), ("JXC", "jxc"), ("",), ("JXC", 3)):
        schema.task_selectors = invalid
        with pytest.raises(ValueError):
            schema.accepted_task_selectors()


def test_incomplete_explicit_ipd_accepts_missing_required_selector():
    source = (_EXAMPLES / "arpes.in").read_text()
    incomplete = source.replace("TASK     ARPES\n", "TASK\n", 1)
    assert source != incomplete
    definition = InputParameters.definition("ARPES")

    with pytest.raises(pp.ParseBaseException):
        definition.read_from_file(io.StringIO(incomplete))

    with pytest.warns(UserWarning, match="Missing non-optional input item.*TASK"):
        parsed = definition.read_from_file(
            io.StringIO(incomplete), allow_incomplete=True
        )
    assert parsed.definition_name == "ARPES"
    assert not parsed.TASK.TASK.is_set()


def test_incomplete_autodetect_does_not_guess_scf_or_shared_task():
    source = (_EXAMPLES / "arpes.in").read_text()
    incomplete = source.replace("TASK     ARPES\n", "TASK\n", 1)
    with pytest.raises(ValueError, match="Cannot autodetect incomplete"):
        InputParameters.from_file(io.StringIO(incomplete), allow_incomplete=True)


def test_incomplete_explicit_task_retains_identity_checks():
    from ..task_definition import TaskMismatchError

    source = (_EXAMPLES / "arpes.in").read_text()
    incomplete = source.replace("TASK     ARPES\n", "TASK\n", 1)
    registry = InputParameters.input_parameters_registry
    task = InputParameters.task_registry["SPLEED"]

    parsed = registry.read_for_task(
        task, io.StringIO(incomplete), allow_incomplete=True
    )
    assert parsed.task_name == "SPLEED"
    assert parsed.requested_task_name == "SPLEED"
    assert not parsed.TASK.TASK.is_set()

    with pytest.raises(TaskMismatchError):
        registry.read_for_task(
            task, io.StringIO(source), allow_incomplete=True
        )


def test_incomplete_mode_still_rejects_invalid_value():
    definition = InputParameters.definition("ARPES")
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string(
            "TASK ARPES\nSPEC_EL NT=not_an_integer\n", allow_incomplete=True
        )


def test_incomplete_parse_uses_validation_reporting_policy():
    from ...common.warnings import DataValidityError, DataValidityWarning
    import warnings

    definition = InputParameters.definition("ARPES")
    source = (_EXAMPLES / "arpes.in").read_text()
    incomplete = source.replace("TASK     ARPES\n", "TASK\n", 1)

    with pytest.warns(DataValidityError, match=r"Missing non-optional input item.*TASK"):
        parsed = definition.read_from_file(
            io.StringIO(incomplete),
            allow_incomplete=True,
            report_invalid="warn",
        )
    assert not parsed.TASK.TASK.is_set()

    with pytest.raises(DataValidityError, match="TASK"):
        definition.read_from_file(
            io.StringIO(incomplete),
            allow_incomplete=True,
            report_invalid="raise",
        )

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always", DataValidityWarning)
        parsed = definition.read_from_file(
            io.StringIO(incomplete),
            allow_incomplete=True,
            report_invalid="ignore",
        )
    assert not any(issubclass(w.category, DataValidityWarning) for w in recorded)
    assert not parsed.TASK.TASK.is_set()


def test_incomplete_parse_does_not_warn_on_failed_parse():
    from ...common.warnings import DataValidityWarning
    import warnings

    definition = InputParameters.definition("ARPES")
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        with pytest.raises(pp.ParseBaseException):
            definition.read_from_string(
                "TASK ARPES\nSPEC_EL NT=not_an_integer\n",
                allow_incomplete=True,
            )

    assert not any(issubclass(w.category, DataValidityWarning) for w in recorded)
