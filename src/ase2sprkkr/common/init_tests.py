""" Routines and classes used in tests """

import sys
import os
from pathlib import Path
import unittest
import numpy as np
import inspect
import sys
import asyncio
from ase.cell import Cell
from contextlib import contextmanager
def patch_package(package, name):
    """ Set the package name for the tests, to make the relative imports working.

    Usage:

    ..code::

      if __package__:
        from .init_tests import TestCase, patch_package
      else:
         from init_tests import TestCase, patch_package
      __package__, __name__ = patch_package(__package__, __name__)

    """
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

class TestCase(unittest.TestCase):
  """ A testcase class with some usefull assertions and a better numpy
  arrays comparison """

  def assertAsyncEqual(self, a, b):
      return self.assertEqual(a, self.runAsync(b))

  def assertAsyncRaises(self, a, b):
      return self.assertRaises(a, self.runAsync(b))

  def assertAlmostEqual(self, a, b, **kwargs):
      np.testing.assert_almost_equal(a,b, **kwargs)

  @staticmethod
  def runAsync(corr):
      return asyncio.run(corr)

  @classmethod
  @contextmanager
  def almost_equal_precision(cls, **kwargs):
      tmp = cls._almost_equal_precision
      cls._almost_equal_precision = kwargs
      yield
      cls._almost_equal_precision = tmp

  _almost_equal_precision = {}

  def assertDictEqual(self, a, b, msg=''):
      if a.__class__ != b.__class__:
         if msg: msg+='\n'
         raise ValueError(msg + f'Classes {a.__class__} and {b.__class__} are not equal')
      if len(a) != len(b):
         super().assertDictEqual(dict(a), dict(b), msg)
      for (ai, av),(bi, bv) in zip(a.items(), b.items()):
         self.assertEqual(ai, bi, 'Dict keys are not equal')
         self.assertEqual(av, bv, 'Dict values are not equal')

  def setUp(self):
      """ Register numpy array for the equality """

      def testfce(fce, msg='', **kwargs):
        def np_array_equal(a, b, msg=msg, **kwar):
          try:
            fce(a,b, **kwargs, **kwar)
          except AssertionError as e:
            if msg:
               msg = msg + '\n' + str(e)
               raise self.failureException(msg)
            else:
               raise
        return np_array_equal

      assert_equals = testfce(np.testing.assert_equal)
      assert_almost_equals = testfce(np.testing.assert_almost_equal)

      def arr_testfce(a,b,msg, **kwargs):
         """ assert_almost_equal does not work for non-numeric dtypes """
         if a.dtype == 'O':
            return assert_equals(a,b,msg)
         if a.dtype.names:
            for i in range(len(a.dtype)):
                if a.dtype[1] == 'O':
                   return assert_equals(a,b,msg)
         kwargs.update(self._almost_equal_precision)
         return assert_almost_equals(a,b,msg, **kwargs)

      self.addTypeEqualityFunc(
         np.ndarray,
         arr_testfce
      )

      self.addTypeEqualityFunc(
         dict,
         self.assertDictEqual
      )


      self.addTypeEqualityFunc(
         Cell,
         lambda a,b,msg: arr_testfce(np.array(a), np.array(b), msg)
      )
