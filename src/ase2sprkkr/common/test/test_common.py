import asyncio
import inspect
import os
import tempfile
import warnings

import pytest

from ..backward_compatibility import ExceptionGroup, _FallbackExceptionGroup
from ..decorators import add_to_signature
from ..file_utils import FilePath
from ..process_output_reader import AsyncioFileReader
from ..warnings import (
    DataValidityError,
    DataValidityErrors,
    DataValidityWarning,
    ValidationResult,
    catch_warnings_of,
)


class TestCommon:
    def test_validation_result(self):
        warning = DataValidityWarning("warning")
        error = DataValidityError("error")
        assert isinstance(warning, ValidationResult)
        assert isinstance(error, ValidationResult)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ValidationResult.emit([warning, error])
        assert [warning, error] == [item.message for item in caught]

        with pytest.raises(TypeError):
            ValidationResult.emit(UserWarning("not a validation result"))

    def test_validation_collection_does_not_capture_other_warnings(self):
        validation = DataValidityWarning("validation")
        error = DataValidityError("error")

        def emit():
            ValidationResult.emit(validation)
            warnings.warn("ordinary", UserWarning)
            raise error

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            results = ValidationResult.collect(emit)

        assert [validation, error] == results
        assert ["ordinary"] == [str(item.message) for item in caught]

    def test_catch_warnings_of_catches_only_requested_category(self):
        validation = DataValidityWarning("validation")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with catch_warnings_of(ValidationResult) as selected:
                warnings.warn(validation)
                warnings.warn("ordinary", UserWarning)

        assert [validation] == selected
        assert ["ordinary"] == [str(item.message) for item in caught]

    def test_validation_reporting_warns_before_raising_errors(self):
        first = DataValidityError("first")
        later = DataValidityWarning("later")
        second = DataValidityError("second")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(DataValidityErrors) as raised:
                ValidationResult.report([first, later, second], "raise")

        assert isinstance(raised.value, DataValidityError)
        assert isinstance(raised.value, ExceptionGroup)
        assert (first, second) == raised.value.exceptions
        assert isinstance(raised.value.derive((first,)), DataValidityErrors)
        assert [later] == [item.message for item in caught]

        with pytest.raises(DataValidityError) as raised:
            ValidationResult.report([first], "raise")
        assert raised.value is first

    def test_fallback_exception_group_can_split_nested_groups(self):
        first = ValueError("first")
        second = TypeError("second")
        nested = _FallbackExceptionGroup("nested", (second,))
        group = _FallbackExceptionGroup("all", (first, nested))

        selected = group.subgroup(ValueError)
        assert selected.exceptions == (first,)

        selected, remainder = group.split(TypeError)
        assert selected.exceptions[0].exceptions == (second,)
        assert remainder.exceptions == (first,)

    def test_file_path(self, tmp_path):
        directory = os.fspath(tmp_path)
        path = FilePath("results/output.out", directory)
        absolute = os.path.join(directory, "results", "output.out")

        assert str(path) == "results/output.out"
        assert os.fspath(path) == absolute
        assert os.path.abspath(path) == absolute
        assert path == FilePath(absolute, directory)

    def test_common(self):
        def f(a, b=2, **kwargs):
            return a, b, kwargs

        @add_to_signature(f)
        def ff(c=2, *args, **kwargs):
            return c, f(*args, **kwargs)

        assert ff("a", c="c", d="d", e="e") == ("c", ("a", 2, {"d": "d", "e": "e"}))
        assert ff("a", "b", "c", d="d", e="e") == ("c", ("a", "b", {"d": "d", "e": "e"}))
        assert "a" in inspect.signature(ff).parameters
        assert list(inspect.signature(ff).parameters.values())[0].name == "a"

        @add_to_signature(f, prepend=True)
        def ff(c, *args, **kwargs):
            return c, f(*args, **kwargs)

        assert ff("c", "a", d="d", e="e") == ("c", ("a", 2, {"d": "d", "e": "e"}))

        def f(c=1, d=5):
            return c, d

        @add_to_signature(f, prepend=True)
        def ff(a=2, d="d", *args, **kwargs):
            return a, f(d=d, *args, **kwargs)

        assert ff("a", "c") == ("a", (1, "c"))

        @add_to_signature(f, prepend=True)
        def ff(a=2, *args, d="d", **kwargs):
            return a, f(d=d, *args, **kwargs)

        assert ff("a", "c") == ("a", ("c", "d"))

        def f(*args):
            return args

        @add_to_signature(f)
        def ff(a, *args):
            return a, f(*args)

        assert ff("a", "c", "d") == ("a", ("c", "d"))

        @add_to_signature(lambda: None)
        def ff(*args, c, **kwargs):
            return c

        assert ff(c=1) == 1

        @add_to_signature(lambda x: x)
        def ff(*args, c, **kwargs):
            return c, args

        assert ff(2, c=1) == (1, (2,))

        @add_to_signature(lambda x: x)
        def ff(*args, c, **kwargs):
            return c, kwargs

        assert ff(x=2, c=1) == (1, {"x": 2})

        @add_to_signature(lambda: None)
        def ff(a, *args, c, **kwargs):
            return a, c

        assert ff(2, c=1) == (2, 1)

        def f(a=1, b=2):
            return a, b

        @add_to_signature(f, excluding="b")
        def ff(c=1, *args, **kwargs):
            return c, f(b=5, *args, **kwargs)

        assert ff(3, 7) == (7, (3, 5))
        assert ff(3, c=7) == (7, (3, 5))
        assert ff(a=3, c=7) == (7, (3, 5))

    def test_asyncio_file_reader(self):
        async def check(reader):
            assert await reader.readline() == b"First line\n"
            assert await reader.readline() == b"Second line\n"
            assert await reader.readuntil(b"is") == b"Hi, this"
            assert await reader.readuntil(b"is") == b" is"
            assert await reader.readuntil(b"!!!") == b" a very long file!! really!!!"
            assert await reader.readuntil(b"123456") == b" 12345123456"
            with pytest.raises(asyncio.IncompleteReadError):
                await reader.readuntil(b"123456")

        with tempfile.TemporaryFile(mode="w+b") as tfile:
            tfile.write(b"First line\nSecond line\nHi, this is a very long file!! really!!! 12345123456789")
            tfile.seek(0)
            reader = AsyncioFileReader(tfile, buffersize=5)
            asyncio.run(check(reader))

    def test_cached_property(self):
        from ..decorators import cached_property

        class A:
            u = 1

            @cached_property
            def a(self):
                self.__class__.u += 1
                return self.__class__.u

        a = A()
        assert a.a == 2
        assert a.a == 2
        a.a = 5
        assert a.a == 5
        assert a.a == 5
        del a.a
        assert a.a == 3
        assert a.a == 3
        a = A()
        assert a.a == 4
        assert a.a == 4

        class A:
            u = 1

            @cached_property
            def a(self):
                self.__class__.u += 1
                return self.__class__.u

            @a.setter
            def a(self, value):
                self.__dict__["a"] = value + 10

            @a.deleter
            def a(self):
                self.__dict__["a"] = 25

        a = A()
        assert a.a == 2
        assert a.a == 2
        a.a = 5
        assert a.a == 15
        assert a.a == 15
        del a.a
        assert a.a == 25
        assert a.a == 25
        a = A()
        assert a.a == 3
        assert a.a == 3

        def modifier(function):
            def wrap(self):
                return function(self)

            return cached_property(wrap)

        class A:
            u = 1

            @modifier
            def a(self):
                self.__class__.u += 1
                return self.__class__.u

            @modifier
            def b(self):
                self.__class__.u += 1
                return self.__class__.u

        a = A()
        assert a.a == 2
        assert a.a == 2
        assert a.b == 3
        assert a.b == 3
        assert a.a == 2
