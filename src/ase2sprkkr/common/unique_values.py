from collections import namedtuple
import numpy as np

class UniqueValuesMapping:

  def __init__(self, mapping, value_to_class_id=None):
      """
      Parameters
      ----------
      mapping: np.ndarray
        Array of equivalence class members
        members[id] = <eq class id>

      reverse: dict
        Mapping { value: <eq class id> }
      """
      self.mapping = mapping
      self.value_to_class_id = value_to_class_id

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
        means "any class". The will be assigned to the first available class of
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

  noneMapping = namedtuple('NoneMapping', ['merge'])(lambda x: Create(values, reverse_mapping=False))
