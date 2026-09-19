"""Definitions of the supported configuration-item repetition modes."""

from enum import Enum, nonmember
from typing import Callable, TYPE_CHECKING, Union
from pyparsing import DelimitedList
import inspect
import numpy as np

from .parsing_results import ArrayKey, DefArrayKey, DefDictKey, DictKey, IgnoredKey, Key, RepeatedKey
from .grammar import Forward

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


def _exact_number(x):
        return x,x


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


class VariableRepeatedItemGrammar(RepeatedItemGrammar):
    """ Grammar constructor for an item with variable number of items """

    def __init__(self, count):
        self.count = count
        if isinstance(count, str):
            self.names = [ count ]
            self.fn = _exact_number
        elif callable(count):
            sig = inspect.signature(count)
            self.names = list(sig.parameters)
            self.fn = count
        else:
            raise TypeError("repeated_count must be an int, str, or callable.")
        self.args = { i:None for i in self.names }
        self.forward = None

    def call(self, *args):
        out = self.fn(*args)
        if isinstance(out, int):
            return out, out
        return out

    def limits(self, item):
        container = item._container
        return self.fn(*(container[name]() for name in self.names))

    def apply_hooks(self, container):

        for name in self.names[:-1]:

            def parse_action(s, l, t, name=name):
                self.args[name] = t[0][1]

            def grammar_hook(grammar, parse_action=parse_action):
                grammar.add_parse_action(parse_action)

            container[name].add_grammar_hook(grammar_hook)

        def parse_action(s, l, t):
            self.args[self.names[-1]] = t[0][1]
            count = self.call(*self.args.values())
            self.forward << self._delimited_list(self.grammar, self.delimiter, count[0], count[1])

        def grammar_hook(grammar):
            grammar.add_parse_action(parse_action)

        container[self.names[-1]].add_grammar_hook(grammar_hook)

    def grammar(self, grammar, delimiter):

        if not self.forward:
            self.forward = Forward()
        self.grammar = grammar
        self.delimiter = delimiter or ''
        return self.forward

    def copy(self):
        return VariableRepeatedItemGrammar(self.count)
