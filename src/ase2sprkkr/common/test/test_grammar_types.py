import warnings
import pyparsing
import numpy as np
import datetime
from ase.units import Rydberg

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from .. import grammar_types as gt      # NOQA: E402
from ..grammar import generate_grammar  # NOQA: E402


class TestGrammar(TestCase):

  def test_is_the_same_value(self):
      self.assertTrue(gt.Integer.is_the_same_value(1,1))
      self.assertFalse(gt.Integer.is_the_same_value(1,0))
      self.assertTrue(gt.Real.is_the_same_value(1.,1.))
      self.assertFalse(gt.Real.is_the_same_value(1.,0.))
      self.assertTrue(gt.String.is_the_same_value("1111","1111"))
      self.assertFalse(gt.String.is_the_same_value("1112","1111"))
      self.assertTrue(gt.Array.is_the_same_value(np.array([1,2,3]),np.array([1,2,3])))
      self.assertFalse(gt.Array.is_the_same_value(np.array([1,2,3]),np.array([1,2,4])))
      self.assertTrue(gt.Array.is_the_same_value([1,2,3],[1,2,3]))
      self.assertFalse(gt.Array.is_the_same_value([1,2,3],[1,2,4]))
      # TODO this is not working. However, it should not be needed, so i let it as it for now
      # self.assertTrue(gt.Array.is_the_same_value(np.array([1,2,np.array([1,2])]),np.array([1,2,np.array([1,2])])))
      # self.assertFalse(gt.Array.is_the_same_value(np.array([1,2,np.array([1,2])]),np.array([1,2,np.array([1,3])])))
      # self.assertTrue(gt.Array.is_the_same_value(np.array([1,2,(np.array([1,2]),1)]),np.array([1,2,(np.array([1,2]),1)])))
      # self.assertFalse(gt.Array.is_the_same_value(np.array([1,2,(np.array([1,2]),1)]),np.array([1,2,(np.array([1,3]),1)])))

  def test_types(self):
    with generate_grammar():
       self._test_types()

  def _test_types(self):
    Error = object()

    class Invalid:

      def __init__(self, v):
          self.value = v

    def test(val, res):
      # try:
        g = type.grammar()
        try:
          out = g.parseString(str(val), True)
          assert len(out) == 1
          out = out[0]
          self.assertEqual(out, res, "{} should be {} and is {} for type {}".format(val, res, out, type.__class__.__name__))
        except (ValueError, pyparsing.ParseException) as e:
          assert res is Error, f"'{val}' should be validated as '{res}'" \
              f" for type {type.__class__.__name__}, an error {e} have been given"
          test_invalid(val)
          return
        assert type.validate(out)
        out = g.parseString(val, True)
        assert len(out) == 1
        out = out[0]
        self.assertEqual(out, res, "{} should be {} and is {} for type {} after input and output".format(val, res, out, type.__class__.__name__))
      # except Exception as e:
      #  if isinstance(e, AssertionError):
      #    raise
      #  raise AssertionError(f"Error '{val}' should be validated to '{res}', {str(e)} occurs")

    def test_invalid(val):
      try:
        with warnings.catch_warnings():
          warnings.simplefilter("ignore")
          val = type.convert(val)
          self.assertTrue(type.validate(val) is not True)
      except ValueError:
        pass
      else:
        self.assertTrue(False)

    def test_warning(val, res):
        with warnings.catch_warnings(record = True) as warning_list:
            val = type.convert(val)
            self.assertEqual(val, res)
            self.assertTrue(type.validate(val))
        self.assertEqual(True, bool(warning_list))

    def test_valid(val):
      val = type.convert(val)
      self.assertTrue(type.validate(val) is True)

    type = gt.Integer()
    for val, res in [
        ('1', 1),
                    ]:
        test(val, res)

    type = gt.Integer(min=-5, max=15)
    for val, res in [
        (-101, Error),
    ]:
        test(val, res)

    type = gt.Integer()
    for val, res in [
        ('1', 1),
        ('a', Error),
        (1.1, Error),
        ('1.1 a', Error),
        ('-1', -1),
        ('-111', -111),
        ('999', 999),
                    ]:
        test(val, res)
    for v in [1.1, 'aaaa', (1,2,3)]:
        test_invalid(v)
    for v,r in [(1., 1)]:
        test_warning(v,r)

    type = gt.Integer(min=-5, max=15)
    for val, res in [
        ('1', 1),
        ('a', Error),
        (1.1, Error),
        ('1.1 a', Error),
        ('-1', -1),
        ('15', 15),
        (16, Error),
        (-6, Error),
        ('-5', -5)
                    ]:
        test(val, res)
    for v in [-10, 1.5, 'aaaa', (1,2,3)]:
        test_invalid(v)
    for v,r in [(1., 1)]:
        test_warning(v,r)

    type = gt.Real()
    for val, res in [
        ('1', 1.0),
        ('a', Error),
        ('-1.1', -1.1),
        ('1.1 a', Error),
        ('-1e-2', -1e-2)
                    ]:
        test(val, res)
    for v in ['aaaa', (1,2,3)]:
        test_invalid(v)

    type = gt.Real(min=-5., max=10.)
    for val, res in [
        ('10', 10.0),
        ('10.0001', Error),
        ('a', Error),
        ('-5.1', Error),
        ('-5', -5.),
        ('1.1 a', Error),
        ('-1e-2', -1e-2)
                    ]:
        test(val, res)
    for v in [15., 'aaaa', (1,2,3)]:
        test_invalid(v)
    for v,r in [(1, 1.0)]:
        test_warning(v,r)

    type = gt.String()
    for val, res in [
         ('a', 'a'),
         ('a a', Error),
         ('1','1')
                    ]:
         test(val, res)
    for v in [1, 15., (1,2,3)]:
        test_invalid(v)

    type = gt.Energy()
    for val, res in [
         ('1', 1.0),
         ('Ry', Error),
         ('1 Ry',1.0),
         ('1 eV', 1.0 / Rydberg),
                    ]:
         test(val, res)
    for v in ['aaaa', (1,2,3)]:
        test_invalid(v)
    for v,r in [(1, 1.0),(15, 15.0)]:
        test_warning(v,r)

    type = gt.SetOf(int)
    for val, res in [
        ('a{}', Error),
        ('{a}', Error),
        ('{1,2}', np.array([1,2])),
        ('{1,a}', Error)
                    ]:
        test(val, res)
    for v in [[1.5,2,3], 'asasd' ]:
        test_invalid(v)
    for v,r in [([1.,2,3], np.asarray([1,2,3]))]:
        test_warning(v,r)
    for v in [[1,2], np.array([1,2,3,4]) ]:
        test_valid(v)

    type = gt.Array(int, min_length=3)
    for val, res in [
        ('a{} 3 4 4', Error),
        ('{a}', Error),
        ('1 2 3 4', np.array([1,2,3,4])),
        ('1 2', Error),
        ('{1,a}', Error)
                    ]:
        test(val, res)
    for v in [[1.5,2,3], 'asasd', 3 ]:
        test_invalid(v)
    for v,r in [([1.,2,3], np.asarray([1,2,3]))]:
        test_warning(v,r)
    for v in [[1,2,3], np.array([1,2,3,4]) ]:
        test_valid(v)

    type = gt.Array(int, min_length=3, as_list=True)
    for val, res in [
        ('a{} 3 4 4', Error),
        ('{a}', Error),
        ('1 2 3 4', [1,2,3,4]),
        ('1 2', Error),
        ('{1,a}', Error)
                    ]:
        test(val, res)

    type = gt.Sequence(int, str, float)
    for val, res in [
        ('a{}', Error),
        ('{1 a 1.1}', Error),
        ('1 aaa 2.', (1,'aaa', 2.)),
                    ]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float)
    for val, res in [(
        """ X YY ZZZ
        1 dog 2.5
        3 cat 3e-2""", np.array([(1,'dog',2.5), (3,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)])
                    ),
        (""" XX YY ZZZ
        1 dog 2.5
        3 cat 3e-2""",
        Error),
        (""" X YY ZZZ
        1.1 dog 2.5
        3 cat 3e-2""",
        Error),
                    ]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True)
    for val, res in [(
        """ X YY ZZZ
        1 1 dog 2.5
        2 3 cat 3e-2""", np.array([(1,'dog',2.5), (3,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)])
                    )]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, grouping=True)
    for val, res in [(
        """ X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 1 5 dog 55.5
        3 2 6 cat 3e-2""",
        [
          np.array([(1,'dog',2.5), (3,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
          np.array([(1,'dog',2.5)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
          np.array([(5,'dog',55.5), (6,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
        ]
                    ),
        (""" X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        4 1 5 dog 3e-2
        4 2 6 cat 3e-2""", Error),
        (""" X YY ZZZ
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 2 5 dog 3e-2
        3 3 6 cat 3e-2""", Error),
                     ]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, grouping=True, group_size='GSIZE')
    for val, res in [(
        """ X YY ZZZ
        GSIZE       2
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        2 2 1 hat 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 0.001 """,
        np.array([[(1,'dog',2.5), (3,'cat', 3e-2)],
          [(1,'dog',2.5), (1,'hat', 2.5)],
          [(5,'dog',3e-2), (6,'cat', 0.001)]], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
                    ),
        (
        """ X YY ZZZ
        GSIZE       3
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        2 2 1 hat 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 3e-2 """, Error),
        (
        """ X YY ZZZ
        GSIZE       2
        1 1 1 dog 2.5
        1 2 3 cat 3e-2
        2 1 1 dog 2.5
        3 1 5 dog 3e-2
        3 2 6 cat 3e-2 """, Error),
                     ]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float, grouping=True, group_size='GSIZE')
    for val, res in [(
        """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 0.2
        2 6 cat 3e-2 """,
        np.array([[(1,'dog',2.5), (3,'cat', 3e-2)],
          [(1,'dog',2.5), (1,'hat', 2.5)],
          [(5,'dog', 0.2), (6,'cat', 3e-2)]], dtype=[('X', int), ('YY', object), ('ZZZ', float)])
                    ),
        (
        """ X YY ZZZ
        GSIZE       3
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 3e-2
        2 6 cat 3e-2""", Error),
        (
        """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        1 5 dog 3e-2
        2 6 cat 3e-2""", Error),
                    ]:
        test(val, res)

    type = gt.Table(X=int, YY=str, ZZZ=float, grouping=True, group_size='GSIZE', groups_as_list=True)
    for val, res in [(
        """ X YY ZZZ
        GSIZE       2
        1 1 dog 2.5
        2 3 cat 3e-2
        1 1 dog 2.5
        2 1 hat 2.5
        1 5 dog 0.2
        2 6 cat 3e-2 """,
        [
          np.array([(1,'dog',2.5), (3,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
          np.array([(1,'dog',2.5), (1,'hat', 2.5)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
          np.array([(5,'dog', 0.2), (6,'cat', 3e-2)], dtype=[('X', int), ('YY', object), ('ZZZ', float)]),
        ]
                   )]:
        test(val, res)

    type = gt.Table(X=int, numbering='Y')
    for val, res in [(
        """ Y X
        1 1
        2 3""", np.array([[1], [3]])
                    )]:
        test(val, res)

    type = gt.Table(X=int, numbering='Y', flatten=True)
    for val, res in [(
        """ Y X
        1 1
        2 3""", np.array([1, 3])
                    )]:
        test(val, res)

    type = gt.PotMixed.I
    for val, res in [
        ('1', 1),
        ('-2',-2),
        ('-2.4',-2.4),
        ('AHOJ', 'AHOJ'),
        ('T', True),
        ('F', False),
        ('aaaa aaaa', 'aaaa aaaa'),
        ('{1,2,3}', np.array([1,2,3])),
        ('{1.,2.,3.}', np.array([1.,2.,3.]))
                    ]:
        test(val, res)

    type = gt.Mixed.I
    for val, res in [
        ('1', 1),
        ('-2',-2),
        ('-2.4',-2.4),
        ('AHOJ', 'AHOJ'),
        ('', True),
        ('aaaa aaaa', Error),
        ('{1,2,3}', np.array([1,2,3])),
        ('{1.,2.,3.}', np.array([1.,2.,3.]))
                   ]:
        test(val, res)
    for v in [ 'asdasdsad', [1.,3.], 8. ]:
        test_valid(v)

    type = gt.Table(X=int, YY=str, ZZZ=float, numbering=True, format='>20', numbering_format='<4')
    data = """ X YY ZZZ
        1 2 dog 2.5"""
    parsed = type.parse(data)
    data2 = type.string(parsed)
    # propper column widths: numbering + 3 columns + newline, no newline on the end of the string
    self.assertEqual(len(data2), 2 * (4 + 3 * 21 +1) - 1)

    type = gt.Range(float)
    for val, res in [
         ('{40,50}', np.array([40.,50.])),
         ('40', 40.),
         ('{40,50,60}', Error),
         ('{40}', Error),
                    ]:
         test(val, res)
    for v in [[1.,2,3], 'asasd' ]:
        test_invalid(v)
    a = lambda *args: np.asarray(args)
    for v,r in [ (3,3.), ((5,1.),a(5.,1.)), ((7.,3),a(7.,3.)) ]:
        test_warning(v,r)
    for v in [ [1.,3.], 8. ]:
        test_valid(v)

    type = gt.Date()
    for val, res in [
         ('{40,50}', Error),
         ('23.01.2020', datetime.date(2020,1,23) ),
         ('40', Error),
                    ]:
         test(val, res)
