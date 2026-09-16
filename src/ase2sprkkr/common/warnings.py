"""
This module implement simple class that encapsulates the
warning about non-critical errors (or just suspicious data)
in configuration/data files.
"""

import warnings
from collections.abc import Iterable
from contextlib import contextmanager
from typing import Any, Callable, Generator, List, Literal, Optional, Sequence, Type

from .backward_compatibility import ExceptionGroup


ValidationReason = Literal["set", "parse", "save"]
RetainInvalidPolicy = Literal["none", "typed", "all"]
ReportInvalidPolicy = Literal["raise", "warn", "ignore"]


@contextmanager
def catch_warnings_of(
    category: Type[Warning],
) -> Generator[List[Warning], None, None]:
    """Record warnings of ``category`` without swallowing other warnings.

    The standard :func:`warnings.catch_warnings` records every warning which
    passes the active filters; its ``category`` argument only installs another
    filter.  This wrapper separates the requested category and re-emits all
    other recorded warnings after restoring the surrounding warning context.

    The yielded list is populated when the context is left.

    Parameters
    ----------
    category
      Warning base class to collect; warnings outside this class are re-emitted.
    """
    recorded = ()
    selected = []
    try:
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always", category)
            yield selected
    finally:
        for warning in recorded:
            if issubclass(warning.category, category):
                selected.append(warning.message)
            else:
                warnings.warn_explicit(
                    warning.message,
                    warning.category,
                    warning.filename,
                    warning.lineno,
                    source=warning.source,
                )


class ValidationResult(UserWarning):
    """Base class for issues returned by semantic validators."""

    @staticmethod
    def emit(result: Any) -> None:
        """Emit one result, or a flat iterable of validation results.

        ``result`` may be ``None``, one :class:`ValidationResult`, or an
        iterable of them.
        """
        if result is None:
            return

        invalid_result = TypeError(
            "A validator must return None, a ValidationResult, "
            "or an iterable of ValidationResult instances"
        )
        if isinstance(result, ValidationResult):
            results = (result,)
        elif isinstance(result, Iterable) and not isinstance(result, (str, bytes)):
            results = tuple(result)
        else:
            raise invalid_result

        if not all(isinstance(issue, ValidationResult) for issue in results):
            raise invalid_result
        for issue in results:
            warnings.warn(issue, stacklevel=3)

    @staticmethod
    def collect(
        callback: Callable[..., Any], *args: Any, **kwargs: Any
    ) -> List["ValidationResult"]:
        """Run ``callback(*args, **kwargs)`` and return its emitted results."""
        try:
            with catch_warnings_of(ValidationResult) as results:
                callback(*args, **kwargs)
        except ValidationResult as issue:
            results.append(issue)
        return results

    @staticmethod
    def report(
        results: Sequence["ValidationResult"],
        invalid: ReportInvalidPolicy = "raise",
    ) -> None:
        """Report every warning, then handle collected errors as requested.

        With ``invalid='raise'``, reporting is deliberately not interrupted by
        the first :class:`DataValidityError`: all ordinary validation warnings
        are emitted first and the collected errors are raised afterwards.
        For  ``'warn'`` respective ``'ignore'``, both the errors and warnings
        are emitted as warnings, or ignored, respectively.

        ``results`` is the collected validation output and ``invalid`` selects
        the ``"raise"``, ``"warn"`` or ``"ignore"`` reporting policy.
        """
        if invalid == 'ignore':
            return
        if invalid == 'warn':
            for result in results:
                warnings.warn(result, stacklevel=4)
        elif invalid == 'raise':
            errors = []
            for result in results:
                if isinstance(result, DataValidityError):
                    errors.append(result)
                else:
                    warnings.warn(result, stacklevel=4)
            if errors:
                if len(errors) == 1:
                    raise errors[0]
                message = "Validation failed: " + "; ".join(
                    str(error) for error in errors
                )
                raise DataValidityErrors(message, errors)
        else:
            raise ValueError(f"Unknown invalid-value reporting policy: {invalid!r}")

class DataValidityWarning(ValidationResult):
    """
    This Warning should be issued, if there are some invalid data
    or format problems, that should not yield a "hard error" which
    would prevent parsing the data.
    """

    @classmethod
    def warn(cls, out: Any) -> None:
        """Emit a warning whose message is built from ``out``."""
        warnings.warn(cls(out), stacklevel=2)


class DataValidityError(DataValidityWarning):
    """Errors of this class will be considered to be 'Errors'.
    The current action will be interrupted."""

    pass


class DataValidityErrors(DataValidityError, ExceptionGroup):
    """An exception group containing at least one validity error."""

    def __init__(self, message: str, exceptions: Sequence[Exception]) -> None:
        """Initialize the exception-group base explicitly for Python 3.8."""
        ExceptionGroup.__init__(self, message, exceptions)

    def derive(self, exceptions: Sequence[Exception]) -> ExceptionGroup:
        """Preserve this type while the derived group has a validity error."""
        if any(isinstance(error, DataValidityError) for error in exceptions):
            return type(self)(self.message, exceptions)
        return ExceptionGroup(self.message, exceptions)


class InvalidValuePolicy:
    """Retention and reporting policy for invalid proposed values."""

    def __init__(
        self,
        why: ValidationReason,
        retain_invalid: Optional[RetainInvalidPolicy] = None,
        report_invalid: Optional[ReportInvalidPolicy] = None,
    ) -> None:
        """Create policies for validation phase ``why``.

        ``retain_invalid`` is ``"none"``, ``"typed"`` or ``"all"``;
        ``report_invalid`` is ``"raise"``, ``"warn"`` or ``"ignore"``.
        Phase-dependent defaults are used when either argument is ``None``.
        """
        if why not in ("set", "parse", "save"):
            raise ValueError(f"Unknown validation reason: {why!r}")
        if retain_invalid is None:
            retain_invalid = "typed" if why == "parse" else "none"
        if report_invalid is None:
            report_invalid = "warn" if why == "parse" else "raise"
        if retain_invalid not in ("none", "typed", "all"):
            raise ValueError(
                f"Unknown invalid-value retention policy: {retain_invalid!r}"
            )
        if report_invalid not in ("raise", "warn", "ignore"):
            raise ValueError(
                f"Unknown invalid-value reporting policy: {report_invalid!r}"
            )
        self.why = why
        self.retain_invalid = retain_invalid
        self.report_invalid = report_invalid
        self.results = []
        self.rollback = False

    @property
    def retain_typed(self) -> bool:
        return self.retain_invalid in ("typed", "all")

    @property
    def retain_all(self) -> bool:
        return self.retain_invalid == "all"

    def add(self, result: ValidationResult) -> None:
        """Record a retained validation ``result``."""
        self.results.append(result)

    def discard(self, result: ValidationResult) -> None:
        """Record an unstaged invalid value and preserve raising atomicity."""
        self.results.append(result)
        if (
            isinstance(result, DataValidityError)
            and self.report_invalid == "raise"
        ):
            self.rollback = True

    def add_semantic(self, results: Sequence[ValidationResult]) -> None:
        """Record semantic ``results`` and request rollback when required."""
        self.results.extend(results)
        if not self.retain_typed and any(
            isinstance(result, DataValidityError) for result in results
        ):
            self.rollback = True

    def report(self) -> None:
        """Apply the configured reporting policy to all recorded results."""
        ValidationResult.report(self.results, self.report_invalid)
