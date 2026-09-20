"""Definitions of the supported configuration-item repetition modes."""

from enum import Enum, nonmember
from numbers import Integral
from typing import Callable, TYPE_CHECKING, Union
from pyparsing import DelimitedList, Empty
import numpy as np

from .dependencies import DependentValue
from .parsing_results import ArrayKey, DefArrayKey, DefDictKey, DictKey, IgnoredKey, Key, RepeatedKey
from .grammar import Forward
from .warnings import DataValidityError

if TYPE_CHECKING:
    from .options import Option


class Repeated(Enum):
    """Describe how repeated configuration items are parsed and stored."""

    @nonmember
    class Type(Enum):
        """Storage used for a repeated item."""

        NO = None
        ARRAY = np.ndarray
        LIST = list
        DICT = dict

        def __bool__(self):
            return self.value

    @nonmember
    class Numbering(Enum):
        """Numbering used by an item name."""

        NO = 0
        YES = 1
        WITH_DEFAULT = 2

        def __bool__(self):
            return self.value > 0

        @property
        def has_default(self):
            return self.value == 2

    def __init__(
        self,
        type: Type,
        key_type: Union[Key, callable] = Key.NONE,
        is_numbered: Numbering = Numbering.NO,
        has_header: bool = True,
    ):
        self.type = type
        self.key_type = key_type
        self.is_numbered = is_numbered
        self.has_header = has_header

    def __bool__(self):
        return self.type != self.Type.NO

    @property
    def grammar_numbering(self):
        """Numbering syntax accepted by the parser."""
        return self.is_numbered

    @property
    def numbering_condition(self):
        """Condition controlling output numbering, if there is one."""
        return None

    def numbering_for(self, option: "Option") -> "Repeated.Numbering":
        """Return the numbering to use when writing ``option``."""
        return self.is_numbered

    @staticmethod
    def NUMBERED_IF(
        condition: Callable[["Option"], bool]
    ) -> "RepeatedNumberedIf":
        """Create a dense repetition numbered when ``condition(option)`` is true."""
        return RepeatedNumberedIf(condition)

    @classmethod
    def create(cls, value, default):
        if value is False:
            return cls.NO
        if value is True:
            return default
        if isinstance(value, str):
            return cls[value]
        return value

    @property
    def is_array(self):
        return self.type in (self.Type.ARRAY, self.Type.LIST)

    @property
    def is_dict(self):
        return self.type == self.Type.DICT

    NO = Type.NO
    IGNORED = (Type.NO, IgnoredKey)
    REPEATED = (Type.LIST, RepeatedKey)
    ARRAY = (Type.ARRAY, RepeatedKey)
    LIST_SECTION = (Type.LIST, Key.NONE, Numbering.NO, False)
    DICT_SECTION = (Type.DICT, Key.NONE, Numbering.NO, False)
    NUMBERED = (Type.ARRAY, ArrayKey, Numbering.YES)
    DICT = (Type.DICT, DictKey, Numbering.YES)
    DEFAULTDICT = (Type.DICT, DefDictKey, Numbering.WITH_DEFAULT)


class RepeatedNumberedIf:
    """Dense repetition whose output numbering is selected at runtime.

    The parser accepts both ``NAME`` and ``NAME<number>`` because the option
    controlling the condition may occur later in the input. Values are always
    stored as an array; ``condition(option)`` only controls validation
    and the spelling used for output.
    """

    type = Repeated.Type.ARRAY
    key_type = DefArrayKey
    is_numbered = Repeated.Numbering.YES
    grammar_numbering = Repeated.Numbering.WITH_DEFAULT
    has_header = True

    def __init__(self, condition: Callable[["Option"], bool]) -> None:
        """Store the one-argument runtime numbering ``condition``."""
        if not callable(condition):
            raise TypeError("The NUMBERED_IF condition has to be callable")
        self.numbering_condition = condition

    def __bool__(self):
        return True

    @property
    def is_array(self):
        return True

    @property
    def is_dict(self):
        return False

    def numbering_for(self, option: "Option") -> Repeated.Numbering:
        """Evaluate numbering for runtime ``option``."""
        return (
            Repeated.Numbering.YES
            if self.numbering_condition(option)
            else Repeated.Numbering.NO
        )


class RepeatedItemGrammar:
    """ Base class for repetition handling grammar constructors """
    @staticmethod
    def create(repeated, count):
        if not repeated:
            return NotRepeatedItemGrammar.I
        if count is None:
            return AnyRepeatedItemGrammar.I
        if isinstance(count, int):
            return ConstRepeatedItemGrammar(count)
        if isinstance(count, RepeatedItemGrammar):
            return count
        return VariableRepeatedItemGrammar(count)

    def apply_hooks(self, container):
        """ Apply any grammar hooks to the container - use for dynamic parsing. """
        pass

    def grammar(self, grammar, delimiter):
        """ Returns a grammar for repeated item """
        raise NotImplementedError()

    def copy(self):
        return self

    def validate(self, item, count, why):
        """ Validate the number of items in the container """
        pass

    @staticmethod
    def _delimited_list(grammar, delimiter, min=None, max=None):
        if min == max == 0:
            return Empty().set_parse_action(lambda _: [])
        if min == max == 1:
            return grammar.copy().add_parse_action(lambda x: x.as_list())
        return DelimitedList(grammar, delimiter, min=min, max=max).set_parse_action(lambda x: x.as_list())

class NotRepeatedItemGrammar(RepeatedItemGrammar):
    """ Crammar constructor for ordinaryu not-repeated item"""

    def grammar(self, grammar, delimiter):
        return grammar

    def apply_hooks(self, container):
        pass

NotRepeatedItemGrammar.I = NotRepeatedItemGrammar()


class AnyRepeatedItemGrammar(RepeatedItemGrammar):
    """ Grammar constructor for item with any number of repetition """

    def grammar(self, grammar, delimiter):
        return self._delimited_list(grammar, delimiter or '')

AnyRepeatedItemGrammar.I = AnyRepeatedItemGrammar()

class ValidatedRepeatedItemGrammar(RepeatedItemGrammar):

    def limits(self, item):
        raise NotImplementedError()

    def validate(self, item, count, why):
        lower, upper = self.limits(item)

        if lower is not None and count < lower:
            if upper is None:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur at least {lower} times"
                )
            elif lower == upper:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur exactly {lower} times"
                )
            else:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur between "
                    f"{lower} and {upper} times"                )

        if upper is not None and count > upper:
            if lower is None:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur at most {upper} times"
                )
            elif lower == upper:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur exactly {lower} times"
                )
            else:
                return DataValidityError(
                    f"Item {item._get_path()} has to occur between "
                    f"{lower} and {upper} times"
                )


class ConstRepeatedItemGrammar(ValidatedRepeatedItemGrammar):
    """ Grammar constructor for item with const number of items """

    def __init__(self, count):
        self.count = count, count if isinstance(count, int) else count

    def limits(self, item):
        return self.count

    def grammar(self, grammar, delimiter):
        return self._delimited_list(grammar, delimiter or '', self.count[0], self.count[1])


class VariableRepeatedItemGrammar(ValidatedRepeatedItemGrammar):
    """ Grammar constructor for an item with variable number of items """

    def __init__(self, count):
        self.dependency = DependentValue.create(count)
        self.count = self.dependency
        self.forward = None

    @staticmethod
    def _limits(value):
        if isinstance(value, Integral):
            value = int(value)
            return value, value
        return value

    def limits(self, item):
        return self._limits(self.dependency.runtime_value(item))

    def apply_hooks(self, container):
        self.dependency.bind(container, self._set_count)

    def _set_count(self, value):
        lower, upper = self._limits(value)
        self.forward << self._delimited_list(
            self._item_grammar, self.delimiter, lower, upper
        )

    def grammar(self, grammar, delimiter):
        if not self.dependency.hooks:
            missing = ", ".join(self.dependency.paths)
            raise KeyError(f"No repeated_count dependencies found: {missing}")

        if not self.forward:
            self.forward = Forward()
        self._item_grammar = grammar
        self.delimiter = delimiter or ''
        return self.forward

    def copy(self):
        return VariableRepeatedItemGrammar(self.dependency)
