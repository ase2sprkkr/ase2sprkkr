"""Values that explicitly bypass configuration validation."""

from __future__ import annotations

import warnings
from typing import Any, Optional, TextIO, TYPE_CHECKING

from .warnings import DataValidityError

if TYPE_CHECKING:
    from .grammar_types import GrammarType


class DangerousValue:
    """A value that bypasses the normal option validation."""

    def __init__(
        self,
        value: Any,
        value_type: Optional["GrammarType"] = None,
        validate: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        value
          A value to be stored.

        value_type
          A grammar type that the value should satisfy. If it is ``None``,
          the only requirement is that the value can be stringified.

        validate
          Whether to validate the value. Values produced by the grammar have
          already been validated and can pass ``False``.
        """
        if validate:
            if value_type:
                with warnings.catch_warnings():
                    warnings.simplefilter("error", DataValidityError)
                    value = value_type.convert(value)
                    value_type.validate(value)
            else:
                value = str(value)
        self.value = value
        self.value_type = value_type

    def __call__(self) -> Any:
        """Return the actual value."""
        return self.value

    def write_value(self, file: TextIO) -> None:
        """Write the stored value to the text ``file``."""
        if self.value_type:
            self.value_type.write(file, self.value)
        else:
            file.write(str(self.value))
