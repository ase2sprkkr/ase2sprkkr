"""
Classes, that represents various value types that can appear in the configuration and problem definitionfiles.

Each grammar type can both parse string containing a value of a given type, and to create the string containing a given value.
"""
from ..grammar import generate_grammar

context =  generate_grammar()
context.__enter__()
#it ensures that the generated grammar will have the correct whitespaces

#will be initialized later
type_from_type_map = {}

from .grammar_type import *
from .basic import *
from .arrays import *
from .mixed import *
from .data import *


type_from_type_map = {
    float  : Real.I,
    complex: Complex.I,
    int  : Integer.I,
    bool : Bool.I,
    str  : String.I
}
""" The standard grammar_types for python types.

The value type can be given by a standard python type, this map maps the
python type for the appropriate grammar_type class.
"""

type_from_set_map = {
    float: set_of_reals,
    int  : set_of_integers,
}
""" Map the python type of a collection member to a grammar type of the collection.

Only canonical types are expected, see :meth:`ase2sprkkr.common.grammar_types.normalize_type`
"""

recognized_set_types = ( list, tuple, np.ndarray )
""" The types, that are recognized as 'list of values' and so that will
be accepted as values for array_like type (e.g. :class:`Array` or :class:`SetOf`). """

#some cleanup
context.__exit__(None, None, None)
del context
