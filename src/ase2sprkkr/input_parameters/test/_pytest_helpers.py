import numpy as np
import pyparsing as pp
import pytest


def as_array(value):
    return np.atleast_1d(value)


def parse(text, grammar):
    if not isinstance(grammar, pp.ParserElement):
        grammar = grammar()
    return grammar.parse_string(text, True)


def assert_equal(actual, expected):
    np.testing.assert_equal(actual, expected)


def assert_parse(text, expected, grammar):
    result = parse(text, grammar).as_list()
    assert len(result) == 1
    actual = result[0]
    if hasattr(actual, "to_dict"):
        actual = actual.to_dict()
    assert_equal(actual, expected)


def assert_not_valid(text, grammar):
    with pytest.raises(pp.ParseBaseException):
        parse(text, grammar)
