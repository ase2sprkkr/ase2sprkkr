import re

import pytest

from ..output_printing import PrintAll, PrintMatching, output_listeners
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
    assert info and len(listeners) == 1 and isinstance(listeners[0], PrintAll)


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
