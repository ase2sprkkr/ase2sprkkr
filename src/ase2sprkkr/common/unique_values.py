""" UniqueValuesMapping: the class for solving equivalence classes on a collection of objects. """
from __future__ import annotations

from collections import namedtuple, defaultdict
from collections.abc import Iterable
import numpy as np
from typing import Union, Dict


class UniqueValuesMapping:
  """ A class, that can map a collection of (possible non-unique) values to a set
      of unique identifiers. It effectively makes the classes of equivalence
      between indexes of the input array.

      The instances of the class can be merged to distinct the values, that are
      the same according to one criterion, but distinct on the other.

      .. doctest::

        >>> UniqueValuesMapping.from_values([1,4,1]).mapping
        array([1, 2, 1])
        >>> UniqueValuesMapping.from_values([int, int, str]).mapping
        array([1, 1, 2])
        >>> UniqueValuesMapping.from_values([1,4,1]).value_to_class_id
        {1: 1, 4: 2}
        >>> UniqueValuesMapping.from_values([1,4,1,1]).merge([1,1,2,1]).mapping
        array([1, 2, 3, 1])
  """

  def __init__(self, mapping:List, value_to_class_id:Dict=None):
      """
      Parameters
      ----------
      mapping: Union[np.ndarray, list]
        Array of equivalence class members
        members[id] = <eq class id>

      reverse: dict
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
        >>> UniqueValuesMapping([1,4,1]).indexes(start_from = 1)
        {1: [1, 3], 2: [2]}
      """

      indexes = {}
      for i,ec in enumerate(self.mapping):
          indexes.setdefault(ec, []).append(i+start_from)
      return indexes

  def iter_unique(self):
      return self.value_to_class_id.keys()

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
        np.array([0,1,0])
        >>> UniqueValuesMapping.from_values([1.,4.,1.]).value_to_class_id
        {1.:0, 4.:1}
      """
      mapping, reverse = UniqueValuesMapping._create_mapping(values, length)
      return UniqueValuesMapping(mapping, reverse)

  @staticmethod
  def _create_mapping(values, length=None):
      mapping = np.empty(length or len(values), int)
      reverse = {}

      for i,v in enumerate(values):
          if v in reverse:
             tag = reverse[v]
          else:
             tag = len(reverse) + 1
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

  def normalize(self, strict:bool=False):
      """ Replace the names of equivalent classes by the integers.

      Parameters
      ----------
      strict
         If True, the resulting integer names will be from range 1..n, where n is the number
         of equivalence classes.
         If False and the names are already integers in a numpy array, do nothing.

      >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalize().mapping
      np.array([0,1,0])
      >>> UniqueValuesMapping.from_values([(0,2),(0,3),(0,2)]).normalize().value_to_class_id[1]
      (0,3)
      """
      if not strict and isinstance(self.mapping, np.ndarray) and np.issubdtype(np.integer, self.mapping.dtype):
           return
      mapping, reverse = self._create_mapping(self.mapping)
      if self.value_to_class_id is not None:
         self.value_to_class_id = { k: reverse[v] for k,v in self.value_to_class_id.items() }
      self.mapping = mapping
