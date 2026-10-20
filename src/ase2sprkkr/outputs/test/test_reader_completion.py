"""Finalize task output exactly once, preserving error precedence."""

from types import SimpleNamespace

import pytest

from ..task_result import KkrOutputReader


class _Parser:
    def __init__(self, print_output=False):
        self.print_output = print_output
        self.calls = []
        self.raise_on_nonzero = False

    def result(self, output, error, return_code):
        self.calls.append((output, error, return_code))
        if self.raise_on_nonzero and return_code:
            raise ValueError("process failed")
        return output, error, return_code


class _Reader(KkrOutputReader):
    parser_class = _Parser


class _Result:
    def __init__(self):
        self.completions = []

    def complete(self, error, return_code):
        self.completions.append((error, return_code))


def _reader():
    return _Reader(None, None, None)


def test_nonzero_exit_finalizes_only_once_even_if_parser_accepts_it():
    reader = _reader()
    result = _Result()
    assert reader.complete_result(result, "stdout", "stderr", 7) is result
    assert reader.parser.calls == [("stdout", "stderr", 7)]
    assert result.completions == [("stderr", 7)]


def test_successful_run_with_parser_exception_does_not_finalize():
    reader = _reader()
    result = _Result()
    failure = RuntimeError("bad stdout")

    with pytest.raises(RuntimeError) as exc:
        reader.complete_result(result, failure, "stderr", 0)

    assert exc.value is failure
    assert reader.parser.calls == []
    assert result.completions == []


def test_process_failure_outweighs_parser_exception():
    reader = _reader()
    reader.parser.raise_on_nonzero = True
    result = _Result()
    parse_error = RuntimeError("incomplete stream")

    with pytest.raises(ValueError, match="process failed"):
        reader.complete_result(result, parse_error, "stderr", 2)

    assert reader.parser.calls == [(None, "stderr", 2)]
    assert result.completions == []


def test_nonraising_error_parser_preserves_stream_exception():
    reader = _reader()
    result = _Result()
    parse_error = RuntimeError("incomplete stream")

    with pytest.raises(RuntimeError) as exc:
        reader.complete_result(result, parse_error, "stderr", 2)

    assert exc.value is parse_error
    assert reader.parser.calls == [(None, "stderr", 2)]
    assert result.completions == []
