import datetime
import warnings

import numpy as np
import pyparsing
import unyt

from .. import grammar_types as gt
from ..grammar import generate_grammar


def _assert_equal(actual, expected, message=None):
    try:
        np.testing.assert_equal(actual, expected)
    except AssertionError as exc:
        if message:
            raise AssertionError(f"{message}\n{exc}") from exc
        raise


class TestGrammar:
    def test_is_the_same_value(self):
        assert gt.Integer.is_the_same_value(1, 1)
        assert not gt.Integer.is_the_same_value(1, 0)
        assert gt.Real.is_the_same_value(1.0, 1.0)
        assert not gt.Real.is_the_same_value(1.0, 0.0)
        assert gt.String.is_the_same_value("1111", "1111")
        assert not gt.String.is_the_same_value("1112", "1111")
        assert gt.Array.is_the_same_value(np.array([1, 2, 3]), np.array([1, 2, 3]))
        assert not gt.Array.is_the_same_value(np.array([1, 2, 3]), np.array([1, 2, 4]))
        assert gt.Array.is_the_same_value([1, 2, 3], [1, 2, 3])
        assert not gt.Array.is_the_same_value([1, 2, 3], [1, 2, 4])

    def test_types(self):
        with generate_grammar():
            self._test_types()

    def _test_types(self):
        Error = object()
        ParseError = object()

        def test(val, res):
            grammar = type.grammar()
            try:
                out = grammar.parse_string(str(val), True)
                assert len(out) == 1
                out = out[0]
                _assert_equal(
                    out,
                    res,
                    "{} should be {} and is {} for type {}".format(
                        val, res, out, type.__class__.__name__
                    ),
                )
            except (ValueError, pyparsing.ParseException) as exc:
                assert res is Error or res is ParseError, (
                    f"'{val}' should be validated as '{res}'"
                    f" for type {type.__class__.__name__}, an error {exc} have been given"
                )
                if res is Error:
                    test_invalid(val)
                else:
                    assert type.validate(val)
                return
            assert type.validate(out)
            out = grammar.parse_string(val, True)
            assert len(out) == 1
            out = out[0]
            _assert_equal(
                out,
                res,
                "{} should be {} and is {} for type {} after input and output".format(
                    val, res, out, type.__class__.__name__
                ),
            )

        def test_invalid(val):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    converted = type.convert(val)
                    assert type.validate(converted) is not True
            except ValueError:
                return
            raise AssertionError(f"{val!r} should be invalid for {type.__class__.__name__}")

        def test_warning(val, res):
            with warnings.catch_warnings(record=True) as warning_list:
                converted = type.convert(val)
                _assert_equal(converted, res)
                assert type.validate(converted)
            assert warning_list

        def test_valid(val):
            converted = type.convert(val)
            assert type.validate(converted) is True

        type = gt.Integer()
        for val, res in [("1", 1)]:
            test(val, res)

        type = gt.Integer(min=-5, max=15)
        for val, res in [(-101, Error)]:
            test(val, res)

        type = gt.QString()
        for val, res in [("1", "1"), ('"aaa bbb"', "aaa bbb"), ("aaa bbb", ParseError)]:
            test(val, res)

        type = gt.Integer()
        for val, res in [
            ("1", 1),
            ("a", Error),
            (1.1, Error),
            ("1.1 a", Error),
            ("-1", -1),
            ("-111", -111),
            ("999", 999),
        ]:
            test(val, res)
        for v in [1.1, "aaaa", (1, 2, 3)]:
            test_invalid(v)
        for v, r in [(1.0, 1)]:
            test_warning(v, r)

        type = gt.Integer(min=-5, max=15)
        for val, res in [
            ("1", 1),
            ("a", Error),
            (1.1, Error),
            ("1.1 a", Error),
            ("-1", -1),
            ("15", 15),
            (16, Error),
            (-6, Error),
            ("-5", -5),
        ]:
            test(val, res)
        for v in [-10, 1.5, "aaaa", (1, 2, 3)]:
            test_invalid(v)
        for v, r in [(1.0, 1)]:
            test_warning(v, r)

        type = gt.Real()
        for val, res in [("1", 1.0), ("a", Error), ("-1.1", -1.1), ("1.1 a", Error), ("-1e-2", -1e-2)]:
            test(val, res)
        for v in ["aaaa", (1, 2, 3)]:
            test_invalid(v)

        type = gt.Real(min=-5.0, max=10.0)
        for val, res in [
            ("10", 10.0),
            ("10.0001", Error),
            ("a", Error),
            ("-5.1", Error),
            ("-5", -5.0),
            ("1.1 a", Error),
            ("-1e-2", -1e-2),
        ]:
            test(val, res)
        for v in [15.0, "aaaa", (1, 2, 3)]:
            test_invalid(v)
        for v, r in [(1, 1.0)]:
            test_warning(v, r)

        type = gt.String()
        for val, res in [("a", "a"), ("a a", Error), ("1", "1")]:
            test(val, res)
        for v in [1, 15.0, (1, 2, 3)]:
            test_invalid(v)

        type = gt.Energy()
        for val, res in [("1", 1.0 * unyt.Ry), ("Ry", Error), ("1 Ry", 1.0 * unyt.Ry), ("1 eV", 1.0 * unyt.eV)]:
            test(val, res)
        for v in ["aaaa", (1, 2, 3)]:
            test_invalid(v)
        for v, r in [(1, 1.0 * unyt.Ry), (15, 15.0 * unyt.Ry)]:
            test_warning(v, r)

        type = gt.SetOf(int)
        for val, res in [("a{}", Error), ("{a}", Error), ("{1,2}", np.array([1, 2])), ("{1,a}", Error)]:
            test(val, res)
        for val, res in [("{1,2)", Error), ("{1,2", Error), ("(1,2)", np.array([1, 2])), ("[1,2]", np.array([1, 2]))]:
            test(val, res)
        for val, res in [("1,2)", Error), ("{1,2]", Error), ("[1,2]]", Error), ("1,2", np.array([1, 2]))]:
            test(val, res)

        for v in [[1.5, 2, 3], "asasd"]:
            test_invalid(v)
        for v, r in [([1.0, 2, 3], np.asarray([1, 2, 3]))]:
            test_warning(v, r)
        for v in [[1, 2], np.array([1, 2, 3, 4])]:
            test_valid(v)

        type = gt.Array(int, min_length=3)
        for val, res in [
            ("a{} 3 4 4", Error),
            ("{a}", Error),
            ("1 2 3 4", np.array([1, 2, 3, 4])),
            ("1 2", Error),
            ("{1,a}", Error),
        ]:
            test(val, res)
        for v in [[1.5, 2, 3], "asasd", 3]:
            test_invalid(v)
        for v, r in [([1.0, 2, 3], np.asarray([1, 2, 3]))]:
            test_warning(v, r)
        for v in [[1, 2, 3], np.array([1, 2, 3, 4])]:
            test_valid(v)

        type = gt.Array(int, min_length=3, as_list=True)
        for val, res in [
            ("a{} 3 4 4", Error),
            ("{a}", Error),
            ("1 2 3 4", [1, 2, 3, 4]),
            ("1 2", Error),
            ("{1,a}", Error),
        ]:
            test(val, res)

        type = gt.Sequence(int, str, float)
        for val, res in [("a{}", Error), ("{1 a 1.1}", Error), ("1 aaa 2.", (1, "aaa", 2.0))]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float)
        for val, res in [
            (
                """ X YY ZZZ
        1 dog 2.5
        3 cat 3e-2""",
                np.array([(1, "dog", 2.5), (3, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
            ),
            (
                """ XX YY ZZZ
        1 dog 2.5
        3 cat 3e-2""",
                Error,
            ),
            (
                """ X YY ZZZ
        1.1 dog 2.5
        3 cat 3e-2""",
                Error,
            ),
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, repeat_sparse=["ZZZ"])
        for val, res in [
            (
                """ X YY ZZZ
        1 dog 2.5
              3e-2
        2 cat 2.5
              0.8""",
                np.array(
                    [(1, "dog", 2.5), (1, "dog", 3e-2), (2, "cat", 2.5), (2, "cat", 0.8)],
                    dtype=[("X", int), ("YY", object), ("ZZZ", float)],
                ),
            )
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True)
        for val, res in [
            (
                """ X YY ZZZ
        1 1 dog 2.5
        2 3 cat 3e-2""",
                np.array([(1, "dog", 2.5), (3, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
            )
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, grouping=True)
        for val, res in [
            (
                """ X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 1 5 dog 55.5
        3 2 6 cat 3e-2""",
                [
                    np.array([(1, "dog", 2.5), (3, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                    np.array([(1, "dog", 2.5)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                    np.array([(5, "dog", 55.5), (6, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                ],
            ),
            (
                """ X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        4 1 5 dog 3e-2
        4 2 6 cat 3e-2""",
                Error,
            ),
            (
                """ X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 2 5 dog 3e-2
        3 3 6 cat 3e-2""",
                Error,
            ),
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, grouping=True, group_size="GSIZE")
        for val, res in [
            (
                """ X YY ZZZ
        GSIZE       2
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        2 2 1 hat 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 0.001 """,
                np.array(
                    [
                        [(1, "dog", 2.5), (3, "cat", 3e-2)],
                        [(1, "dog", 2.5), (1, "hat", 2.5)],
                        [(5, "dog", 3e-2), (6, "cat", 0.001)],
                    ],
                    dtype=[("X", int), ("YY", object), ("ZZZ", float)],
                ),
            ),
            (
                """ X YY ZZZ
        GSIZE       3
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        2 2 1 hat 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 3e-2 """,
                Error,
            ),
            (
                """ X YY ZZZ
        GSIZE       2
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 3e-2 """,
                Error,
            ),
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, grouping=True, group_size="GSIZE")
        for val, res in [
            (
                """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 0.2
        2 6 cat 3e-2 """,
                np.array(
                    [
                        [(1, "dog", 2.5), (3, "cat", 3e-2)],
                        [(1, "dog", 2.5), (1, "hat", 2.5)],
                        [(5, "dog", 0.2), (6, "cat", 3e-2)],
                    ],
                    dtype=[("X", int), ("YY", object), ("ZZZ", float)],
                ),
            ),
            (
                """ X YY ZZZ
        GSIZE       3
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 3e-2
        2 6 cat 3e-2""",
                Error,
            ),
            (
                """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        1 5 dog 3e-2
        2 6 cat 3e-2""",
                Error,
            ),
        ]:
            test(val, res)

        type = gt.Table(X=int, YY=str, ZZZ=float, grouping=True, group_size="GSIZE", groups_as_list=True)
        for val, res in [
            (
                """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 0.2
        2 6 cat 3e-2 """,
                [
                    np.array([(1, "dog", 2.5), (3, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                    np.array([(1, "dog", 2.5), (1, "hat", 2.5)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                    np.array([(5, "dog", 0.2), (6, "cat", 3e-2)], dtype=[("X", int), ("YY", object), ("ZZZ", float)]),
                ],
            )
        ]:
            test(val, res)

        type = gt.Table(X=int, numbering="Y")
        for val, res in [
            (
                """ Y X
        1 1
        2 3""",
                np.array([[1], [3]]),
            )
        ]:
            test(val, res)

        type = gt.Table(X=int, numbering="Y", flatten=True)
        for val, res in [
            (
                """ Y X
        1 1
        2 3""",
                np.array([1, 3]),
            )
        ]:
            test(val, res)

        type = gt.PotMixed.I
        for val, res in [
            ("1", 1),
            ("-2", -2),
            ("-2.4", -2.4),
            ("AHOJ", "AHOJ"),
            ("T", True),
            ("F", False),
            ("aaaa aaaa", "aaaa aaaa"),
            ("{1,2,3}", np.array([1, 2, 3])),
            ("{1.,2.,3.}", np.array([1.0, 2.0, 3.0])),
        ]:
            test(val, res)

        type = gt.Mixed.I
        for val, res in [
            ("1", 1),
            ("-2", -2),
            ("-2.4", -2.4),
            ("AHOJ", "AHOJ"),
            ("", True),
            ("aaaa aaaa", ParseError),
            ("{1,2,3}", np.array([1, 2, 3])),
            ("{1.,2.,3.}", np.array([1.0, 2.0, 3.0])),
        ]:
            test(val, res)
        for v in ["asdasdsad", [1.0, 3.0], 8.0]:
            test_valid(v)

        type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, format=">20", numbering_format="<4")
        data = """ X YY ZZZ
        1 2 dog 2.5"""
        parsed = type.parse(data)
        data2 = type.string(parsed)
        assert len(data2) == 2 * (4 + 3 * 21 + 1) - 1

        type = gt.Range(float)
        for val, res in [("{40,50}", np.array([40.0, 50.0])), ("40", 40.0), ("{40,50,60}", Error), ("{40}", Error)]:
            test(val, res)
        for v in [[1.0, 2, 3], "asasd"]:
            test_invalid(v)
        array = lambda *args: np.asarray(args)  # noqa: E731
        for v, r in [(3, 3.0), ((5, 1.0), array(5.0, 1.0)), ((7.0, 3), array(7.0, 3.0))]:
            test_warning(v, r)
        for v in [[1.0, 3.0], 8.0]:
            test_valid(v)

        type = gt.Date()
        for val, res in [("{40,50}", Error), ("23.01.2020", datetime.date(2020, 1, 23)), ("40", Error)]:
            test(val, res)
