import sys
import os
from pathlib import Path
import unittest
import numpy as np
import inspect
import sys
import asyncio

def patch_package(package, name):
    """ Set the package name for the test """
    if package and package.count('.') >= 3:
       return package, name
    frame=inspect.stack()[1]
    path = frame.filename
    file = Path(path).resolve()
    current = file.parents[0]
    while file.name != 'ase2sprkkr':
      file = file.parent
    top = str(file.parent)
    sys.path.append(top)
    package=str(current)[len(top)+1:].replace('/','.')
    return package, package + '.' + name.rsplit('.', 1)[-1]

__package__, __name__ = patch_package(__package__,__name__)
from .misc import OrderedDict

class TestCase(unittest.TestCase):

  def assertAsyncEqual(self, a,b):
      return self.assertEqual(a, self.runAsync(b))

  def assertAsyncRaises(self, a,b):
      return self.assertRaises(a, self.runAsync(b))

  @staticmethod
  def runAsync(corr):
      loop = asyncio.get_event_loop()
      task = loop.create_task(corr)
      return loop.run_until_complete(task)

  def setUp(self):
      """ Register numpy array for the equality """

      def testfce(fce, msg=''):
        def np_array_equal(a, b, msg):
          try:
            fce(a,b)
          except AssertionError as e:
            if msg:
               msg = msg + '\n' + str(e)
               raise self.failureException(msg)
            else:
               raise
        return np_array_equal

      assert_equals = testfce(np.testing.assert_equal)
      assert_almost_equals = testfce(np.testing.assert_almost_equal)

      def arr_testfce(a,b,msg):
         """ assert_almost_equal does not work for non-numeric dtypes """
         if a.dtype == 'O':
            return assert_equals(a,b,msg)
         if a.dtype.names:
            for i in range(len(a.dtype)):
                if a.dtype[1] == 'O':
                   return assert_equals(a,b,msg)
         return assert_almost_equals(a,b,msg)

      self.addTypeEqualityFunc(
         np.ndarray,
         arr_testfce
      )

      self.addTypeEqualityFunc(
         OrderedDict,
         assert_equals
      )