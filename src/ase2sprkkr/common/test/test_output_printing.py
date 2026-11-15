import re

import pytest

from ..output_printing import (
    PrintAll, PrintMatching, output_listeners, normalize_output_filter,
    serialize_output_regexes, parse_output_regexes, output_filter_from_gui,
)
from ..process_output_reader import ProcessOutputParser
from ...outputs.readers.scf import ScfOutputParser
from ...outputs.readers.scf import ScfOutputReader
from ...configuration import config


def test_modes():
    assert output_listeners(False) == ((), False)
    assert output_listeners(None) == ((), False)
    assert output_listeners("info") == ((), True)
    listeners, info = output_listeners(True)
    assert not info and len(listeners) == 1 and isinstance(listeners[0], PrintAll)
    listeners, info = output_listeners(r"ERROR")
    assert not info and len(listeners) == 1 and isinstance(listeners[0], PrintMatching)
    listeners, info = output_listeners(("info", r"ERROR", re.compile("WARN")))
    assert info and len(listeners) == 1 and isinstance(listeners[0], PrintMatching)
    assert len(listeners[0].patterns) == 2
    listeners, info = output_listeners((True, "ERROR", "info"))
    assert not info and len(listeners) == 1 and isinstance(listeners[0], PrintAll)


def test_split_lines(capsys):
    (listener,), print_info = output_listeners(r"ERR|WARNING")
    assert print_info is False
    listener("ignore\nER", "out")
    listener("R\nWARN", "out")
    listener("ING", "out")
    listener("ERR in stderr\n", "err")
    listener.flush()
    listener.flush()
    assert capsys.readouterr() == ("ERR\nWARNING", "ERR in stderr\n")


def test_print_all(capsys):
    (listener,), print_info = output_listeners(True)
    assert print_info is False
    listener("hello", "out")
    listener("problem", "err")
    assert capsys.readouterr() == ("hello", "problem")


def test_parser_receives_only_boolean():
    assert ProcessOutputParser(print_info=True).print_info is True
    assert ScfOutputParser(print_info=True).print_info is True
    assert ScfOutputReader(None, None, None, print_info=True).parser.print_info is True


def test_validation():
    with pytest.raises(TypeError):
        output_listeners(22)
    with pytest.raises(TypeError):
        output_listeners(re.compile(b"bad"))
    with pytest.raises(re.error):
        output_listeners("[invalid")


def test_config_accepts_regex_tuple():
    original = config.running.print_output()
    try:
        config.running.print_output.set(r"WARN|ERROR")
        assert config.running.print_output() == r"WARN|ERROR"
        config.running.print_output.set(("info", r"ERROR"))
        assert tuple(config.running.print_output()) == ("info", r"ERROR")
    finally:
        config.running.print_output.set(original)



@pytest.mark.parametrize("raw,mode,regex_count,all_output,info", [
    (None, "none", 0, False, False),
    ("info", "info_regex", 0, False, True),
    (("info", "ERROR", "WARNING"), "info_regex", 2, False, True),
    (("ERROR", "ERROR"), "regex", 1, False, False),
    ((False, "WARNING"), "regex", 1, False, False),
    ((True, "info", "ERROR"), "all", 0, True, False),
    (re.compile("info"), "regex", 1, False, False),
])
def test_normalized_filter_modes(raw, mode, regex_count, all_output, info):
    normalized = normalize_output_filter(raw)
    assert normalized.mode == mode
    assert len(normalized.regexes) == regex_count
    assert normalized.all_output is all_output
    assert normalized.info is info


def test_matcher_uses_search_not_match_and_returns_boolean():
    single = normalize_output_filter("ERROR")
    assert single.matcher("prefix ERROR suffix") is True
    assert single.matcher("ok") is False
    assert type(single.matcher("ERROR")) is bool

    combined = normalize_output_filter(("ERROR", "WARNING"))
    assert combined.matcher("prefix WARNING suffix") is True
    assert combined.matcher("prefix ERROR suffix") is True
    assert combined.matcher("not relevant") is False
    assert normalize_output_filter(None).matcher("anything") is False
    assert normalize_output_filter(True).matcher("anything") is True


def test_overlapping_regexes_print_line_only_once(capsys):
    listener = PrintMatching((re.compile("ERR"), re.compile("ERROR")))
    listener("prefix ERROR suffix\n", "out")
    listener.flush()
    assert capsys.readouterr().out == "prefix ERROR suffix\n"


def test_printmatching_uses_normalized_tuple_predicate(capsys):
    listener = PrintMatching((re.compile("ERROR"), re.compile("WARNING")))
    listener("hello\nERROR\nwarning\nWARNING\n", "out")
    assert capsys.readouterr().out == "ERROR\nWARNING\n"


@pytest.mark.parametrize("sources", [
    (r"ERROR|WARNING",),
    (r"ERROR", r"WARNING"),
    (r"ERROR", r'say "hello"', r"\d+"),
    (r'"quoted',),
    (r"",),
    (r'foo"|bar', r'"WARNING'),
])
def test_regex_text_serialization_is_unambiguous(sources):
    patterns = tuple(re.compile(source) for source in sources)
    text = serialize_output_regexes(patterns)
    parsed = parse_output_regexes(text)
    # The backslash escape for a lone leading quote changes the source
    # spelling but not the matching behavior.
    assert len(parsed) == len(patterns)
    for expected, restored in zip(patterns, parsed):
        for sample in ('"quoted words', "ERROR WARNING", 'say "hello"', "foo", "bar", ""):
            assert bool(expected.search(sample)) == bool(restored.search(sample))


def test_gui_filter_tuple_roundtrip_and_special_info():
    assert output_filter_from_gui("info_regex", "") == "info"
    assert output_filter_from_gui("regex", r"ERROR|WARNING") == r"ERROR|WARNING"
    assert output_filter_from_gui("regex", '"ERROR"|"WARNING"') == ("ERROR", "WARNING")
    assert output_filter_from_gui("info_regex", '"ERROR"|"WARNING"') == (
        "info", "ERROR", "WARNING"
    )
    assert output_filter_from_gui("regex", "info") == "(?:info)"
    value = ("info", re.compile("ERROR", re.IGNORECASE), re.compile(r"\d+"))
    normalized = normalize_output_filter(value)
    saved = output_filter_from_gui(normalized.mode, normalized.regex_text)
    restored = normalize_output_filter(saved)
    assert restored.info
    assert restored.matcher("an error")
    assert restored.matcher("number 42")


def test_output_regex_tuple_quote_escaping():
    parsed = parse_output_regexes('"a""b"|"c|d"')
    assert tuple(p.pattern for p in parsed) == ('a"b', 'c|d')
    assert serialize_output_regexes(parsed) == '"a""b"|"c|d"'
    with pytest.raises(ValueError):
        parse_output_regexes('"unterminated')
    with pytest.raises(ValueError):
        parse_output_regexes('"a"|')
    with pytest.raises(ValueError):
        parse_output_regexes('"a"oops')
    with pytest.raises(re.error):
        parse_output_regexes('"["')


def test_invalid_filter_value_and_binary_pattern():
    with pytest.raises(TypeError):
        normalize_output_filter(re.compile(b"bad"))
    with pytest.raises(TypeError):
        normalize_output_filter(22)
    with pytest.raises(ValueError, match="Invalid output regular expression"):
        output_filter_from_gui("regex", "[")
    with pytest.raises(ValueError):
        output_filter_from_gui("regex", "")
