""" Routines and classes used in tests """

import sys
from pathlib import Path
import numpy as np
import inspect
import asyncio
from ase.cell import Cell
from contextlib import contextmanager
import pytest
import tempfile
import unittest
import os


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
    try:
        import ase2sprkkr # NOQA
    except ImportError:
        sys.path.append(top)
    package=str(current)[len(top) + 1:].replace('/','.')
    return package, package + '.' + name.rsplit('.', 1)[-1]


__package__, __name__ = patch_package(__package__,__name__)


class extdict(dict):

    def __call__(self, **kwargs):
        out = self.copy()
        out.update(kwargs)
        return out


class TestCase:
  """ A testcase class with some usefull assertions and a better numpy
  arrays comparison """

  _print_output = '-v' in sys.argv or '--verbose' in sys.argv
  _calc_args = extdict(
       directory = False, input_file = 'output_test_calc.inp',  # empty_spheres=False,
       output_file = 'output_test_calc.out', potential_file ='output_test_calc.pot', print_output=_print_output,
       mpi = 'auto', options = {'NKTAB': 5, 'NE': 20},
       empty_spheres = False
  )

  @pytest.fixture
  def temporary_dir(self):
      with tempfile.TemporaryDirectory() as d:
          TestCase._calc_args['directory'] = d
          self.dirname = d
          yield
          TestCase._calc_args['directory'] = False
          del self.dirname

  @classmethod
  def calc_args(cls, **kwargs):
      if not kwargs:
          return cls._calc_args
      out = cls._calc_args.copy()
      if 'options' in kwargs:
          kwargs['options'].update(cls._calc_args['options'])
      out.update(kwargs)
      return out

  def run_sprkkr(self):
     return os.environ.get('DO_NOT_RUN_SPRKKR', '') == ''

  def assertAsyncEqual(self, a, b):
      return self.assertEqual(a, self.runAsync(b))

  def assertAsyncRaises(self, a, b):
      with pytest.raises(a):
          self.runAsync(b)

  def assertRaises(self, a, b=None):
      if b is None:
          return pytest.raises(a)
      with pytest.raises(a):
          b()

  def assertAlmostEqual(self, a, b, **kwargs):
      np.testing.assert_almost_equal(a,b, **kwargs)

  def assertIsNone(self, a):
      assert a is None

  @staticmethod
  def runAsync(corr):
      return asyncio.run(corr)

  def assertTrue(self, val):
      assert val

  def assertFalse(self, val):
      assert not val

  def assertEqual(self, a, b, msg=None):
      assertion.assertEqual(a,b, msg=msg)

  def assertNotEqual(self, a, b):
      with pytest.raises(AssertionError):
          assertion.assertEqual(a,b)

  @classmethod
  @contextmanager
  def almost_equal_precision(cls, **kwargs):
      tmp = assertion._almost_equal_precision
      assertion._almost_equal_precision = kwargs
      yield
      assertion._almost_equal_precision = tmp


assertion = unittest.TestCase('__init__')
assertion._almost_equal_precision = {}


def testfce(fce, msg='', **kwargs):
    def np_array_equal(a, b, msg=msg, **kwar):
        try:
          fce(a,b, **kwargs, **kwar)
        except AssertionError as e:
          if msg:
             msg = msg + '\n' + str(e)
             raise assertion.failureException(msg)
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
    kwargs.update(assertion._almost_equal_precision)
    return assert_almost_equals(a,b,msg, **kwargs)


assertion.addTypeEqualityFunc(
    np.ndarray,
    arr_testfce
)


def assertDictEqual(a, b, msg=''):
    if a.__class__ != b.__class__:
       if msg:
           msg+='\n'
       raise ValueError(msg + f'Classes {a.__class__} and {b.__class__} are not equal')
    if len(a) != len(b):
       super().assertDictEqual(dict(a), dict(b), msg)
    for (ai, av),(bi, bv) in zip(a.items(), b.items()):
       assertion.assertEqual(ai, bi, 'Dict keys are not equal')
       assertion.assertEqual(av, bv, 'Dict values are not equal')


def assertListEqual(a, b, msg=''):

   def message(error):
       return msg + ': ' + error if msg else error

   assert a.__class__ is b.__class__, message('A list is expected')
   assert len(a) == len(b), message('The lists should have the same lengths, they have '
                                f'the lengths {len(a)} and {len(b)} respectivelly')
   for i, vals in enumerate(zip(a,b)):
       try:
          assertion.assertEqual(vals[0], vals[1])
       except Exception as e:
          raise AssertionError(message(f'The {i}th value of the lists differs: {e}'))


assertion.addTypeEqualityFunc(list, assertListEqual)
assertion.addTypeEqualityFunc(dict, assertDictEqual)
assertion.addTypeEqualityFunc(Cell,
    lambda a,b,msg: arr_testfce(np.array(a), np.array(b), msg)
)
