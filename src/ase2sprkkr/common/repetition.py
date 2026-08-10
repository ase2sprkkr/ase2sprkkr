"""Definitions of the supported configuration-item repetition modes."""

from enum import Enum, nonmember
from typing import Union

import numpy as np

from .parsing_results import ArrayKey, DefArrayKey, DefDictKey, DictKey, IgnoredKey, Key, RepeatedKey


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

    def numbering_for(self, option):
        """Return the numbering to use when writing ``option``."""
        return self.is_numbered

    @staticmethod
    def NUMBERED_IF(condition):
        """Create a dense repeated value numbered when ``condition(option)`` is true."""
        return RepeatedNumberedIf(condition)

    @classmethod
    def create(cls, value, grammar_type=None):
        if value is False:
            return cls.NO
        if value is True:
            if hasattr(grammar_type, "is_numpy_array") and grammar_type.is_numpy_array:
                return cls.ARRAY
            return cls.REPEATED
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
    stored as an array; ``condition(option)`` only controls validation and the
    spelling used for output.
    """

    type = Repeated.Type.ARRAY
    key_type = DefArrayKey
    is_numbered = Repeated.Numbering.YES
    grammar_numbering = Repeated.Numbering.WITH_DEFAULT
    has_header = True

    def __init__(self, condition):
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

    def numbering_for(self, option):
        return Repeated.Numbering.YES if self.numbering_condition(option) else Repeated.Numbering.NO
