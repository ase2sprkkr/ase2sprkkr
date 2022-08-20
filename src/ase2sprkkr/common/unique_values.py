""" UniqueValuesMapping: the class for solving equivalence classes on a collection of objects. """

from collections import namedtuple, defaultdict
import numpy as np

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

  def __init__(self, mapping, value_to_class_id=None):
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
        >>> UniqueValuesMapping.from_values([1,4,1]).indexes(start_from = 1)
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

  @staticmethod
  def from_values(values, length=None, none_tuples=False):
      """
      values: iterable
        Values to find the equivalence classes

      length: int
        Length of values - provide it, if len(values) is not available

      none_tuples: bool
        If True, the items of the values are 2-items tuples. Then, None in the second one
        means "any class". Such one will be assigned to the first available class of
        the first tuple's value.
      """

      mapping = np.empty(length or len(values), int)
      reverse = {}

      for i,v in enumerate(values):
          if v in reverse:
             tag = reverse[v]
          elif none_tuples and (v[0], None) in reverse:
             tag = reverse[(v[0], None)]
             if v[1] is not None:
                del reverse[(v[0], None)]
                reverse[v] = tag
          else:
             tag = len(reverse) + 1
             reverse[v] = tag
          mapping[i] = tag

      return UniqueValuesMapping(mapping, reverse)

  def merge(self, other):
      """ Merge two sets. Resulting UniqueValues uses integers as keys"""
      return self.from_values(zip(self.mapping, other), length = len(self.mapping), none_tuples=True)
