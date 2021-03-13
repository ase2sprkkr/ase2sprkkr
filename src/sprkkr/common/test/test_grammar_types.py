if __package__:
   from .init_tests import TestCase
else:
   from init_tests import TestCase, patch_package
   __package__ = patch_package()

import pyparsing
from .. import grammar_types as gt
import numpy as np
from ase.units import Rydberg

class GrammarTest(TestCase):

  def test_types(self):
    Error = "ERROR"

    type = gt.Integer()



    def test(val, res):
      #try:
        try:
          g = type.grammar()
          out = g.parseString(val, True)
          assert len(out) == 1
          out = out[0]
          self.assertEqual(out, res, "{} should be {} and is {} for type {}".format(val, res, out, type.__class__.__name__))
        except (ValueError, pyparsing.ParseException) as e:
          assert res is Error, "{} should be validated as {} for type {}".format(val, res, type.__class__.__name__)
          return
        assert type.validate(out)
        out = g.parseString(val, True)
        assert len(out) == 1
        out = out[0]
        self.assertEqual(out, res, "{} should be {} and is {} for type {} after input and output".format(val, res, out, type.__class__.__name__))
      #except Exception as e:
      #  if isinstance(e, AssertionError):
      #    raise
      #  raise AssertionError(f"Error '{val}' should be validated to '{res}', {str(e)} occurs")


    for val, res in [
        ('1', 1),
        ('a', Error),
        ('1.1', Error),
        ('1.1 a', Error),
        ('-1', -1)
        ]:
        test(val, res)

    type = gt.Real()
    for val, res in [
        ('1', 1.0),
        ('a', Error),
        ('1.1', 1.1),
        ('1.1 a', Error),
        ('-1e-2', -1e-2)
        ]:
        test(val, res)

    type = gt.String()
    for val, res in [
         ('a', 'a'),
         ('a a', Error),
         ('1','1')
         ]:
         test(val, res)

    type = gt.Energy()
    for val, res in [
         ('1', 1.0),
         ('Ry', Error),
         ('1 Ry',1.0),
         ('1 eV', 1.0/Rydberg),
         ]:
         test(val, res)



    type = gt.SetOf(int)
    for val, res in [
        ('a{}', Error),
        ('{a}', Error),
        ('{1,2}', np.array([1,2])),
        ('{1,a}', Error)
        ]:
        test(val, res)

    type = gt.Array(int, min_length=3)
    for val, res in [
        ('a{} 3 4 4', Error),
        ('{a}', Error),
        ('1 2 3 4', np.array([1,2,3,4])),
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

    type = gt.Mixed.I
    for val, res in [
        ('1', 1),
        ('-2',-2),
        ('-2.4',-2.4),
        ('AHOJ', 'AHOJ'),
        ('T', True),
        ('F', False),
        ('aaaa aaaa', Error),
        ('{1,2,3}', np.array([1,2,3])),
        ('{1.,2.,3.}', np.array([1.,2.,3.]))
        ]:
        test(val, res)

