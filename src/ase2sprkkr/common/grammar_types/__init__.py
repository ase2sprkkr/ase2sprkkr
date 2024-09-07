"""
Classes, that represents various value types that can appear in the configuration and problem definitionfiles.

Each grammar type can both parse string containing a value of a given type, and to create the string containing a given value.
"""
from ..grammar import generate_grammar
import numpy as np

context = generate_grammar()
context.__enter__()
# it ensures that the generated grammar will have the correct whitespaces

# will be initialized later
type_from_type_map = {}

from .grammar_type import *    # NOQA
from .basic import *           # NOQA
from .arrays import *          # NOQA
from .mixed import *           # NOQA
from .data import *            # NOQA


type_from_type_map = {
    float  : Real.I,
    np.float64 : Real.I,
    complex : Complex.I,
    np.complex128 : Complex.I,
    int  : Integer.I,
    np.int32  : Integer.I,
    bool : Bool.I,
    str  : String.I
}
""" The standard grammar_types for python types.

The value type can be given by a standard python type, this map maps the
python type for the appropriate grammar_type class.
"""

type_from_set_map = {
    float: set_of_reals,
    np.float32: set_of_reals,
    int  : set_of_integers,
    np.int32  : set_of_integers,
}
""" Map the python type of a collection member to a grammar type of the collection.

Only canonical types are expected, see :meth:`ase2sprkkr.common.grammar_types.normalize_type`
"""

recognized_set_types = ( list, tuple, np.ndarray )
""" The types, that are recognized as 'list of values' and so that will
be accepted as values for array_like type (e.g. :class:`Array` or :class:`SetOf`). """

# some cleanup
context.__exit__(None, None, None)
del context
