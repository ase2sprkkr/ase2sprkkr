""" Various classes and routines used thorough the package """

import copy
import numpy as np
import operator


def copy_list(src):
    """ Copy list of objects. Each object is copied just once (so
    the number of unique objects in the list is retained """

    mp = {}

    def cpy(i):
        if i in mp:
           return mp[i]
        mp[i] = copy.copy(i)
        return mp[i]

    return [ cpy(i) for i in src ]


try:
  from numba import njit
except ImportError:
  def njit(fce):
      """ Mock the numba JIT compiler, if it is not available """
      return fce


@njit
def numpy_index(array, item):
    """ Returns index of the first occurence of the item in the array
    If numba is installed, the function is accelerated.
    """
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx
    # If no item was found return None, other return types might be a problem due to
    # numbas type inference.


def as_integer(value):
    """ Interpret the value as integer, or raise (even for float, complex etc.) TypeError.i
    >>> as_integer(4)
    4
    >>> as_integer(np.int64(4))
    4
    >>> as_integer(4.0)
    Traceback (most recent call last):
    ...
    TypeError: 'float' object cannot be interpreted as an integer
    """
    return operator.index(value)


def dict_first_item(self):
    """
    Return the first item of the given dict.
    >>> dict_first_item({'a':4, 'c':2, 'h':8})
    4
    """
    return self[next(iter(self))]


def first_non_none(*args):
    for i in args:
        if i is not None:
            return i
    return None
