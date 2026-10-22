"""Grammar types for variant types (that can accepts values of more types) and the types derived from it (e.g. the Range type is in fact
Variant type that can accepts either number or two-item array)"""

import pyparsing as pp
import copy
import unyt

from .grammar_type import GrammarType, compare_numpy_values, recognized_set_types, type_from_type, type_from_value
from ..decorators import cached_property, add_to_signature

from .basic import Energy, Real, Integer, QString, LineString, Flag, Bool, Boolean
from .arrays import SetOf, set_of_integers, set_of_reals


class BaseMixed(GrammarType):
    """
    A variant type - it can hold "anything".
    """

    type = None
    """ The types, that the value can hold. To be redefined in the descendants. """

    string_type = None
    """ Type of string grammar_type to be used.  To be redefined in the descendants. """

    is_the_same_value = staticmethod(compare_numpy_values)

    def _grammar(self, param_name=False):
        return pp.Or((i.grammar(param_name) for i in self.types))

    def get_type(self, value):
        """Select the grammar type corresponding to a Python value."""
        if isinstance(value, str):
            return self.string_type
        if (isinstance(value, recognized_set_types) and not isinstance(value, unyt.unyt_quantity)
                and len(value) and isinstance(value[0], str)):
            for candidate in self.types:
                if isinstance(candidate, SetOf) and isinstance(candidate.type, self.string_type.__class__):
                    return candidate
        return type_from_value(value)

    def _validate(self, value, why="set"):
        if value is None:
            return True
        type = self.get_type(value)
        if type is value:
            return "Can not determine the type of value {}".format(value)
        return type.validate(value, why)

    def write(self, f, val):
        """Output the value to the stream (in the propper format)."""
        type = self.get_type(val)
        if type is val:
            return super().write(f, val)
        type.write(f, val)

    def grammar_name(self):
        return "<mixed>"

    def _string(self, value):
        if isinstance(value, bool):
            return Bool._string(self, value)
        return self.get_type(value).string(value)

    def convert(self, value):
        if value is None:
            return None
        return self.get_type(value).convert(value)

    def copy_value(self, value):
        return copy.deepcopy(value)


class Range(BaseMixed):
    """A range type - it accepts either one value or range of two values of a given type."""

    @add_to_signature(BaseMixed.__init__, prepend=True)
    def __init__(self, type, *args, **kwargs):
        self._type = type_from_type(type)
        super().__init__(*args, **kwargs)

    @cached_property
    def types(self):
        return [self._type, SetOf(self._type, min_length=2, max_length=2)]

    def get_type(self, value):
        return self.types[1 if isinstance(value, recognized_set_types) else 0]


class Mixed(BaseMixed):
    """Unknown input value; nonempty string sets require explicit brackets."""

    string_type = QString.I
    types = (
        Energy.I, Real.I, Integer.I, set_of_integers, set_of_reals,
        Flag.I, QString.I, SetOf(QString.I, min_length=1, require_brackets=True),
    )
    missing_value = lambda self: (True, True, False)


class PotMixed(BaseMixed):
    """Potential-file unknown value: preserve its line-oriented string grammar."""

    string_type = LineString.I
    types = (Energy.I, Real.I, Integer.I, set_of_integers,
             set_of_reals, Bool.I, LineString.I)


class Variant(BaseMixed):
    """Variant type for configuration values."""

    string_type = QString.I
    types = (
        Energy.I, Real.I, Integer.I, set_of_integers, set_of_reals,
        Boolean.I, QString.I, SetOf(QString.I, min_length=1, require_brackets=True),
    )


class CustomMixed(BaseMixed):
    """Mixed type restricted to an explicit set of grammar types.

    Unlike :class:`BaseMixed`, Python values must be dispatched only to one of
    the configured alternatives. In particular, the string alternative is
    used only for actual strings; it must not silently stringify mappings or
    other unsupported objects.
    """

    @add_to_signature(BaseMixed.__init__, prepend=True)
    def __init__(self, types, string_type=LineString.I, *args, **kwargs):
        self.types = tuple(types)
        self.string_type = string_type
        super().__init__(*args, **kwargs)

    def get_type(self, value):
        if isinstance(value, str):
            if self.string_type not in self.types:
                raise ValueError(f"String values are not allowed for {self}")
            return self.string_type

        if isinstance(value, recognized_set_types):
            for candidate in self.types:
                if not getattr(candidate, "array_access", False):
                    continue
                try:
                    converted = candidate.convert(value)
                    candidate.validate(converted)
                except (TypeError, ValueError, KeyError):
                    continue
                return candidate
            raise ValueError(f"Unsupported sequence value {value!r} for {self}")

        try:
            inferred = type_from_value(value)
        except (TypeError, ValueError, KeyError) as exc:
            raise ValueError(f"Unsupported value {value!r} for {self}") from exc

        for candidate in self.types:
            if candidate is inferred or candidate.__class__ is inferred.__class__:
                return candidate
        raise ValueError(f"Value {value!r} has a type not allowed for {self}")


mixed = Mixed.I = Mixed()  # NOQA: E741
""" A standard grammar type instance for variant (mixed) in input files """
pot_mixed = PotMixed.I = PotMixed()  # NOQA: E741
""" A standard grammar type instance for variant (mixed) values in potential files """
