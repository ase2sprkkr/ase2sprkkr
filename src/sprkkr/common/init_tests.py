import sys
import os
from pathlib import Path
import unittest
import numpy as np
import inspect
import sys

def patch_package():
    """ Set the package name for the test """
    frame=inspect.stack()[1]
    path = frame.filename
    file = Path(path).resolve()
    current = file.parents[0]
    while file.name != 'sprkkr':
      file = file.parent
    top = str(file.parent)
    sys.path.append(top)
    return str(current)[len(top)+1:].replace('/','.')


class TestCase(unittest.TestCase):

  def setUp(self):
      """ Register numpy array for the equality """
      def np_array_equal(a, b, msg='NumPy arrays are not equal'):
          try:
            np.testing.assert_array_equal(a,b)
          except AssertionError:
             raise self.failureException(msg)

      self.addTypeEqualityFunc(
         np.ndarray,
         np_array_equal
      )
