"""Tests for live process-output filtering and SCF summary selection."""

import re

import pytest

from ..output_printing import OutputPrinter, normalize_print_output
from ...outputs.readers.scf import ScfOutputParser
from ...configuration import config


def test_print_output_modes(capsys):
    for mode in (False, None, "info"):
        printer = OutputPrinter(mode)
        assert not printer.enabled
        printer("anything\n", "out")
        printer.flush()
        assert capsys.readouterr() == ("", "")

    printer = OutputPrinter(True)
    assert printer.enabled
    printer("partial", "out")
    printer(" line\n", "out")
    printer("error\n", "err")
    printer.flush()
    captured = capsys.readouterr()
    assert captured.out == "partial line\n"
    assert captured.err == "error\n"


def test_regex_matches_complete_lines_across_chunks(capsys):
    printer = OutputPrinter(r"ERR|WARNING")
    printer("ignored\nER", "out")
    printer("R: failure\nnot a match\nWARN", "out")
    printer("ING without newline", "err")
    printer("ING at end", "out")
    printer.flush()
    printer.flush()
    captured = capsys.readouterr()
    assert captured.out == "ERR: failure\nWARNING at end"
    assert captured.err == ""


def test_regex_filters_stderr_independently(capsys):
    printer = OutputPrinter(re.compile("warning", re.I))
    printer("war", "err")
    printer("ning\n", "err")
    printer("no\n", "out")
    assert capsys.readouterr() == ("", "warning\n")


def test_info_can_be_combined_with_regex(capsys):
    mode = ("info", re.compile("ERR"))
    show_all, info, regexes = normalize_print_output(mode)
    assert not show_all and info and len(regexes) == 1
    assert ScfOutputParser(mode).print_info is True
    assert ScfOutputParser(mode).print_output is True

    printer = OutputPrinter(mode)
    printer("ERR\n", "out")
    printer.flush()
    assert capsys.readouterr().out == "ERR\n"

    info_only = ScfOutputParser(("info",))
    assert info_only.print_info and not info_only.print_output
    assert not ScfOutputParser(False).print_info


def test_invalid_print_output_fails_early():
    with pytest.raises(TypeError, match="print_output"):
        OutputPrinter(23)
    with pytest.raises(TypeError, match="bytes"):
        OutputPrinter(re.compile(b"bad"))
    with pytest.raises(re.error):
        OutputPrinter("[unclosed")


def test_configuration_accepts_regex_and_tuple():
    original = config.running.print_output
    try:
        config.running.print_output = r"ERR|WARNING"
        assert config.running.print_output == r"ERR|WARNING"
        config.running.print_output = ("info", r"ERR|WARNING")
        assert tuple(config.running.print_output) == ("info", r"ERR|WARNING")
    finally:
        config.running.print_output = original
