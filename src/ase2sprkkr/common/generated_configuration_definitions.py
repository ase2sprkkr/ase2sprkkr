"""Definition of value in configuration/output files,
that is generated from other values"""

from .configuration_definitions import RealItemDefinition
from .value_definitions import InheritingValueModifier
from .decorators import add_to_signature
from .options import Option
import copy
import numpy as np
import pyparsing as pp
from .parsing_results import ValidateKey
from typing import Any, Callable, List, Optional, Union
from .warnings import DataValidityError
from .configuration_transaction import ConfigurationTransaction


_NO_KEY = object()


def _validate_length_values(
    definition: RealItemDefinition,
    container: Any,
    length: Union[int, bool] = False,
) -> List[DataValidityError]:
    """Validate source lengths in ``container`` against optional ``length``."""
    first = None
    issues = []
    for name in definition._length_of:
        source = container[name]
        if source.is_dangerous():
            continue
        value = source()
        current_length = None if value is None else len(value)
        if length is False:
            length = current_length
            first = name
        elif length != current_length:
            if first is None:
                message = f"Lengths of {name} should be {length}"
            else:
                message = f"Lengths of {first} and {name} should not differ"
            issues.append(DataValidityError(message))
    return issues


def _validate_length(
    option: Option, container: Any, why: str
) -> List[DataValidityError]:
    """Validate the generated length ``option`` for phase ``why``."""
    length = False
    if why == "parse":
        parsed = getattr(container, "_parsed_values", None)
        if parsed is not None:
            length = getattr(parsed, "checks", {}).get(
                option.name, False
            )
    return _validate_length_values(option._definition, container, length)


class BaseGeneratedValueDefinition(RealItemDefinition):
    """Base class for all generated values. It just set
    that it is generated."""

    is_generated = True
    """ Generated value - the value is computed from other values """
    is_stored = False
    """ This property sets, that this Value/Option is generated. """
    is_validated = False
    """ By default, there is nothing to check on generated values """

    result_class = Option
    """ The generated Values creates :class:`Option` """

    _grammar = None
    """ Do not generated grammar, since this item is not readed, but computed from other values. """

    item_type = "generated value"

    def __repr__(self):
        return f"<{self.name} (generated)>"

    def _generic_info(self):
        return f"Calculated value {self.name}"

    def data_description(self, verbose: Union[bool, str, int] = False, show_hidden: bool = False, prefix: str = ""):
        return ""


class GeneratedValueDefinition(BaseGeneratedValueDefinition):
    """Definition of a value computed from its runtime section.

    ``getter`` receives ``(section)`` for ordinary access and
    ``(section, key)`` for indexed access. ``setter``, when supplied, receives
    ``(section, value)`` or ``(section, value, key)``. A setter that needs the
    active transaction can obtain it using
    :meth:`ConfigurationTransaction.current`.
    """

    @add_to_signature(RealItemDefinition.__init__, prepend=True)
    def __init__(
        self,
        name: str,
        getter: Callable[..., Any],
        setter: Optional[Callable[..., None]] = None,
        **kwargs: Any,
    ) -> None:
        """Define ``name`` using its runtime ``getter`` and optional ``setter``."""
        super().__init__(name, **kwargs)
        self.getter = getter
        self._setter = setter

    @property
    def setter(self):
        if not self._setter:
            raise ValueError(f"Setting the value(s) of {self.name} is not allowed")
        return self._setter


class Length(InheritingValueModifier):
    """Sometimes, the length of some array should appear in the config file"""

    is_generated = True
    is_validated = True
    is_stored = True

    type = int

    _NO_DEFAULT = object()

    validators = (_validate_length,)

    def __init__(self, *of, default_values=_NO_DEFAULT):
        self._length_of = of
        self._has_default_values = default_values is not self._NO_DEFAULT
        if self._has_default_values:
            if not isinstance(default_values, (list, tuple)) or (len(of) == 1 and len(default_values) != 1):
                default_values = (default_values,) * len(of)
            else:
                assert len(default_values) == len(of)
        self._default_values = default_values if self._has_default_values else None

    def modify_definition(self, what):
        what._length_of = self._length_of
        what._length_has_default_values = self._has_default_values
        what._length_defaults = self._default_values
        return super().modify_definition(what)

    def copy(self, **kwargs: Any) -> "Length":
        """Copy this modifier while preserving length metadata overrides."""
        out = super().copy(**kwargs)
        out._length_of = self._length_of
        out._length_has_default_values = self._length_has_default_values
        out._length_defaults = self._length_defaults
        return out

    def stage_generated(
        self,
        option: Option,
        transaction: ConfigurationTransaction,
        value: Any,
        key: Any = _NO_KEY,
    ) -> None:
        """Resize source arrays to ``value`` within ``transaction``.

        ``key`` is accepted for the generated-value protocol but length
        assignment always applies to whole values.
        """
        if not self._length_has_default_values:
            raise DataValidityError(
                f"{self.get_path()} is the length of {' and '.join(self._length_of)}: so it is generated automatically "
                "and can not be set."
            )
        for i, d in zip(self._length_of, self._length_defaults):
            source = option._container[i]
            if value is None or value == 0:
                source.stage(transaction, None)
                continue
            val = source()
            if val is None:
                source.stage(transaction, [d] * value)
                continue
            ln = len(val)
            if ln > value:
                source.stage(transaction, val[:value])
            else:
                if isinstance(val, np.ndarray):
                    extra = len(val.shape) - 1
                    shape = ((0, value - ln),) + ((0, 0),) * extra
                    val = np.pad(val, shape)
                    val[ln:] = d
                else:
                    val = val + [d] * (value - ln)
                source.stage(transaction, val)

    def getter(self, container: Any, key: Any = None) -> Optional[int]:
        """Return the length of the first source value in ``container``."""
        val = container[self._length_of[0]]()
        if val is None:
            return None
        return len(val)

    def validate_parsed(self, option: Option) -> None:
        """Check whether ``option`` was present as required by parse conditions."""
        section = option._container
        parsed = getattr(section, "_parsed_values", None)
        if parsed is None:
            return
        present = self.name in parsed or self.name in getattr(parsed, "checks", ())
        allowed = self.allowed(section)
        if present and not allowed:
            raise pp.ParseException(f"Option {self.get_path()} is not allowed for the current configuration")
        if not present and allowed and not self.is_optional:
            raise pp.ParseException(f"Required option {self.get_path()} is missing")

    def _create_grammar(self, allow_dangerous):
        """Add check for the length"""

        def mark_for_validation(tokens: Any) -> tuple:
            """Mark parsed length ``tokens`` for deferred validation."""
            return ValidateKey(tokens[0]), tokens[1]

        out = super()._create_grammar()
        return out.set_parse_action(mark_for_validation)


class NumpyViewDefinition(BaseGeneratedValueDefinition):
    """
    Values described by this description are possibly reshaped views into a large
    "raw data array"

    Parameters
    ----------
    name
      Name of the resulting variable

    data
      The source variable, from which the data are taken

    selector
      The selector, which select the data to be viewed. Any slice
      or simple numpy index is allowed.

    shape
      The data will be reshaped to given shape. The dimension can be given
      either by number, or by names of other container variables. E.g.
      ``('NE', 5)``

    transpose
     Transpose the source data before returning or indexing.
     If reorder is given, this settings has no effect.

    reorder
     Reorder the axes after reshaping. Argument is the array
     of axes order, e.g. (2,0,1) shifts the last axis to be the first.

    transform_key
     Transform function for the keys idnexing the array
     (e.g. the string name can be transformed to a propper numerical index)

    plot
     PlotInfo object that defines how the results are plotted
    """

    def data_description(self, verbose: Union[bool, str] = False, show_hidden=False, prefix: str = ""):
        if self.shape:
            shape = f"({self.shape})"
        else:
            shape = ""

        out = f"{prefix}{self.name} : view of {self.data}{shape}"
        return out

    @add_to_signature(RealItemDefinition.__init__, prepend=True)
    def __init__(
        self,
        name,
        data,
        selector=slice(None),
        shape=None,
        transpose=False,
        reorder=None,
        transform_key=None,
        *args,
        **kwargs,
    ):
        super().__init__(name, *args, **kwargs)
        self.selector = selector
        self.shape = shape
        self.data = data
        self.transform_key = transform_key
        if reorder:
            self.reorder = reorder
        else:
            self.reorder = transpose

    def determine_shape(self, container: Any) -> tuple:
        """Return the shape of the resulting array, possibly computed using
        properties of the other values in ``container``.

        """

        def get(i):
            if isinstance(i, str):
                return container[i]()
            return i

        return tuple([get(i) for i in self.shape])

    def source(self, container: Any) -> Any:
        """Return the selected source view in ``container``."""
        data = container[self.data]()
        if callable(self.selector):
            out = self.selector(data, container)
        else:
            out = data[self.selector]
        if self.shape:
            out.shape = self.determine_shape(container)
        if self.reorder:
            out = np.transpose(out, axes=self.reorder if self.reorder is not True else None)
        return out

    def getter(self, container: Any, key: Any = None) -> Any:
        """Return the complete view or the item selected by ``key``."""
        out = self.source(container)
        if key is not None:
            if self.transform_key:
                key = self.transform_key(key, container)
            out = out[key]
        return out

    def setter(
        self,
        container: Any,
        value: Any,
        key: Any = slice(None),
    ) -> None:
        """Assign ``value`` at ``key`` in the mutable source view.

        Mutable view writes are deliberately immediate and cannot be rolled
        back.
        """
        if self.transform_key:
            key = self.transform_key(key, container)
        self.source(container)[key] = value

    def copy_value(self, value, all_values=False):
        return copy.copy(value)
