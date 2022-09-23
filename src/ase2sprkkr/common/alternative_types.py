""" This module contains helper functions and data to deal with the problem, that
numpy arrays contains their own datatypes and thus e.g. an ``int`` may not be what
it seems """

import numpy as np
import datetime as datetime

types_alternatives = {
  int: (np.int64,),
  float: (np.float64,),
  complex: (np.complex128,),
  bool: (np.bool_,),
  datetime.datetime: (np.datetime64, ),
}
"Map of types and sets of alternative types"""

normalize_type_map = {}
""" Mapping of alternative types to the 'canonical ones'. """

#fill the map
for i in types_alternatives:
  for j in types_alternatives[i]:
      normalize_type_map[j] = i
del i,j

allowed_types = {
  i : (*j, i) for i,j in types_alternatives.items()
}
""" All types, that are allowed for a given type. I.e., the content of types_alternatives and the primary
type itself """

def normalize_type(type):
    """ Return the 'canonical type' for a given type.

    I.e. it maps numpy internal types to standard python ones.

    doctest:
    >>> normalize_type(np.int64)
    <class 'int'>
    """
    return normalize_type_map.get(type, type)
