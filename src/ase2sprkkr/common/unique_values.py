""" UniqueValuesMapping: the class for solving equivalence classes on a collection of objects. """
from __future__ import annotations

from collections.abc import Iterable
import numpy as np
from typing import Union, Dict, List, Optional


class UniqueValuesMapping:
  """ A class, that can map a collection of (possible non-unique) values to a set
      of unique identifiers. It effectively makes the classes of equivalence
      between indexes of the input array.

      The instances of the class can be merged to distinct the values, that are
      the same according to one criterion, but distinct on the other.

      .. doctest::

        >>> UniqueValuesMapping.from_values([1,4,1]).mapping
        array([1, 2, 1], dtype=int32)
        >>> UniqueValuesMapping.from_values([int, int, str]).mapping
        array([1, 1, 2], dtype=int32)
        >>> UniqueValuesMapping.from_values([1,4,1]).value_to_class_id
        {1: 1, 4: 2}
        >>> UniqueValuesMapping.from_values([1,4,1,1]).merge([1,1,2,1]).mapping
        array([1, 2, 3, 1], dtype=int32)
  """

  def __repr__(self):
      if np.issubdtype(self.mapping.dtype, np.integer):
          v = self.normalized(dtype=False)[0]
      else:
          v = self.mapping
      return f"<UniqueValuesMapping: {v}>"

  def __init__(self, mapping:List, value_to_class_id:Dict=None):
      """
      Parameters
      ----------
      mapping: Union[np.ndarray, list]
        Array of equivalence class members
        members[id] = <eq class id>

      value_to_class_id: dict
        Mapping { value: <eq class id> }
      """

      #: Map from <object index> to <object equivalence class id>.
      self.mapping = mapping
      #: Map from <object> to <object equivalence class id>.
      #: If two mappings are merged, this attribute is not available.
      self.value_to_class_id = value_to_class_id

  def indexes(self, start_from:int=0):
      """
      Returns the dictionary that maps equivalence class id to the list
      of class members indexes.

      Parameters
      ----------
      start_from:
        The indexes are by default zero-based, however they can start with
        the given number (typically with 1).

      ..doctest::
        >>> UniqueValuesMapping([1,4,1]).indexes()
        {1: [0, 2], 4: [1]}
        >>> UniqueValuesMapping([1,4,1]).indexes(start_from = 1)
        {1: [1, 3], 4: [2]}
      """

      indexes = {}
      for i,ec in enumerate(self.mapping):
          indexes.setdefault(ec, []).append(i + start_from)
      return indexes

  def unique_indexes(self):
      """
      Returns the dictionary that maps equivalence class id to the list
      of class members indexes.

      ..doctest::
        >>> UniqueValuesMapping([1,1,4]).unique_indexes()
        [0, 2]
      """
      out = []
      done = set()
      for i, cid in enumerate(self.mapping):
          if cid not in done:
              done.add(cid)
              out.append(i)
      return out

  def iter_unique(self):
      return self.value_to_class_id.keys()

  def unique_items(self):
      return self.value_to_class_id.items()

  def len_of_unique(self):
      return len(self.value_to_class_id)

  def __len__(self):
      return len(self.mapping)

  def __iter__(self):
      return iter(self.mapping)

  @staticmethod
  def from_values(values, length:Optional[int]=None):
      """
      Create equivalence-classes mapping. Unlike the constructor,
      this method tags the values by integers and also compute the reverse
      (value to equivalence class) mapping.

      values: iterable
        Values to find the equivalence classes

      length: int
        Length of values - provide it, if len(values) is not available

      .. doctest::
        >>> UniqueValuesMapping.from_values([1.,4.,1.]).mapping
        array([1, 2, 1], dtype=int32)
        >>> UniqueValuesMapping.from_values([1.,4.,1.]).value_to_class_id
        {1.0: 1, 4.0: 2}
      """
      mapping, reverse = UniqueValuesMapping._create_mapping(values, length)
      return UniqueValuesMapping(mapping, reverse)

  @staticmethod
  def _create_mapping(values, length=None, start_from=1, dtype=np.int32):
      """
      Returns
      -------
      mapping : np.ndarray
            maps the value indexes to equivalence class id

      reverse : dict
            maps equivalence classes to value indexes

      .. doctest::
        >>> UniqueValuesMapping._create_mapping([1.,4.,1.])
        (array([1, 2, 1], dtype=int32), {1.0: 1, 4.0: 2})
      """
      mapping = np.empty(length or len(values), dtype=dtype)
      reverse = {}

      for i,v in enumerate(values):
          if v in reverse:
             tag = reverse[v]
          else:
             tag = len(reverse) + start_from
             reverse[v] = tag
          mapping[i] = tag
      return mapping, reverse

  def merge(self, other):
      """ Merge two sets. Resulting UniqueValues uses integers as keys"""
      return self.from_values(zip(self.mapping, other), length = len(self.mapping))

  def is_equivalent_to(self, mapping:Union[UniqueValuesMapping,Iterable]) -> bool:
      """
      Return, whether the mapping is equal to given another mapping, regardless the actual "names" of the equivalence classes.

      Parameters
      ----------
      mapping
        The other mapping can be given either by instance of this class, or just by any iterable (that returns equivalence class names for the items)

      .. doctest::

        >>> UniqueValuesMapping([1,4,1]).is_equivalent_to([0,1,0])
        True
        >>> UniqueValuesMapping([1,4,1]).is_equivalent_to([0,0,0])
        False
        >>> UniqueValuesMapping([1,4,1]).is_equivalent_to([0,1,1])
        False
        >>> UniqueValuesMapping([1,4,1]).is_equivalent_to([5,3,5])
        True
        >>> UniqueValuesMapping([1,4,1]).is_equivalent_to(UniqueValuesMapping.from_values([2,5,2]))
        True
      """
      return self.are_equivalent(self, mapping)

  @staticmethod
  def are_equivalent(a:Union[UniqueValuesMapping,Iterable],b:Union[UniqueValuesMapping,Iterable]) -> bool:
      """
      Return, whether the two mappings are equal, regardless the actual "names" of the equivalence classes.

      See :meth:`is_equivalent<ase2sprkkr.common.unique_values_mapping.UniqueValuesMapping.is_equivalent_to>`
      """

      mp = {}
      js = set()
      for i,j in zip(a,b):
          if i in mp:
             if mp[i] != j:
                return False
          else:
             if j in js:
                return False
             js.add(j)
             mp[i]=j
      return True

  def normalized(self, start_from=1, strict:bool=True, dtype=None):
      """ Map the class ids to integers

      Parameters
      ----------
      strict
         If True, the resulting integer names will be from range (start_from)..(n+start_from-1),
         where n is the number of equivalence classes.
         If False and the names are already integers in a numpy array, do nothing.

      start_from
         Number the equivalent classes starting from.

      Returns
      -------
      mapping : np.ndarray
         Array of integer starting from start_from, denotes the equivalence classes for the values,
         It holds, that ``mappind[index] == equivalence_class``
      reverse : dict
         Dict ``{ equivalence_class : value }``

      .. doctest::

        >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalized()
        (array([1, 2, 1], dtype=int32), {1: 1, 2: 2})
        >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalized(start_from=0)
        (array([0, 1, 0], dtype=int32), {1: 0, 2: 1})
      """

      if not strict and isinstance(self.mapping, np.ndarray):
          ttype = np.integer if dtype is None else dtype
          if np.issubdtype(ttype, self.mapping.dtype):
              return
      if dtype is False:
          dtype = np.integer
      elif dtype is None:
          dtype = np.int32
      mapping, reverse = self._create_mapping(self.mapping, start_from=start_from, dtype=dtype)
      return mapping, reverse

  def normalize(self, start_from=1, strict:bool=False, dtype=None):
      """ Replace the names of equivalent classes by the integers.

      Parameters
      ----------
      strict
         If True, the resulting integer names will be from range (start_from)..(n+start_from-1),
         where n is the number of equivalence classes.
         If False and the names are already integers in a numpy array, do nothing.

      start_from
         Number the equivalent classes starting from.

      dtype
         dtype of the normalized values. None means ``numpy.int32``, however if not strict,
         any integer type will be sufficient.

      Returns
      -------
      unique_values_mapping
         Return self.

      .. doctest::

        >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalize().mapping
        array([1, 2, 1], dtype=int32)
        >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalize().value_to_class_id[(0,3)]
        2
        >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalize(start_from=0).mapping
        array([0, 1, 0], dtype=int32)
      """
      self.mapping, self.reverse = self.normalized(start_from, strict, dtype)

      if self.value_to_class_id is not None:
         self.value_to_class_id = { k: self.reverse[v] for k,v in self.value_to_class_id.items() }
      return self
