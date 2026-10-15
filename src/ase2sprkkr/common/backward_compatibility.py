"""Backward-compatibility helpers for supported Python versions."""

import enum
from typing import Callable, Optional, Sequence, Tuple, Type, Union

from . import decorators


ExceptionCondition = Union[
    Type[Exception],
    Tuple[Type[Exception], ...],
    Callable[[Exception], bool],
]


class _FallbackExceptionGroup(Exception):
    """Python <3.11 implementation of the used ``ExceptionGroup`` API."""

    def __init__(self, message: str, exceptions: Sequence[Exception]) -> None:
        """Create a group named ``message`` from non-empty ``exceptions``."""
        if not isinstance(message, str):
            raise TypeError("ExceptionGroup message must be a string")
        exceptions = tuple(exceptions)
        if not exceptions:
            raise ValueError("ExceptionGroup must contain at least one exception")
        if not all(isinstance(error, Exception) for error in exceptions):
            raise TypeError("ExceptionGroup can contain only Exception instances")
        super().__init__(message)
        self.message = message
        self.exceptions = exceptions

    def derive(self, exceptions: Sequence[Exception]) -> "_FallbackExceptionGroup":
        """Return the same group type containing ``exceptions``."""
        return type(self)(self.message, exceptions)

    @staticmethod
    def _matches(condition: ExceptionCondition, exception: Exception) -> bool:
        """Return whether ``exception`` satisfies a split condition."""
        if isinstance(condition, type) or (
            isinstance(condition, tuple)
            and all(isinstance(item, type) for item in condition)
        ):
            return isinstance(exception, condition)
        return condition(exception)

    def subgroup(
        self, condition: ExceptionCondition
    ) -> Optional["_FallbackExceptionGroup"]:
        """Return the subgroup matching ``condition``, preserving nesting."""
        return self.split(condition)[0]

    def split(
        self, condition: ExceptionCondition
    ) -> Tuple[
        Optional["_FallbackExceptionGroup"],
        Optional["_FallbackExceptionGroup"],
    ]:
        """Partition this group according to ``condition``."""
        if self._matches(condition, self):
            return self, None

        matching = []
        remaining = []
        for exception in self.exceptions:
            if isinstance(exception, _FallbackExceptionGroup):
                selected, rejected = exception.split(condition)
                if selected is not None:
                    matching.append(selected)
                if rejected is not None:
                    remaining.append(rejected)
            elif self._matches(condition, exception):
                matching.append(exception)
            else:
                remaining.append(exception)

        selected = self.derive(matching) if matching else None
        rejected = self.derive(remaining) if remaining else None
        return selected, rejected

    def __str__(self) -> str:
        count = len(self.exceptions)
        suffix = "exception" if count == 1 else "exceptions"
        return f"{self.message} ({count} sub-{suffix})"


try:
    from builtins import ExceptionGroup
except ImportError:
    ExceptionGroup = _FallbackExceptionGroup


if not hasattr(enum, "nonmember"):

    def _nonmember(cls):
        class Nonmember(decorators.cached_class_property):
            def __getattr__(self, name):
                return getattr(cls, name)

        return Nonmember(lambda: cls)

    enum.nonmember = _nonmember
