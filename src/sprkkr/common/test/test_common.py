if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..misc import add_to_signature
import inspect

class CommonTest(TestCase):

  def test_common(self):

    def f(a, b=2, **kwargs):
        return a,b,kwargs

    @add_to_signature(f)
    def ff(*args, c=2, **kwargs):
        return c,f(*args,**kwargs)

    self.assertEqual( ('c',('a',2, {'d':'d', 'e':'e'})),
                       ff('a',c='c', d='d', e='e'))

    self.assertTrue('a' in inspect.signature(ff).parameters)
    self.assertEqual('a', list(inspect.signature(ff).parameters.values())[0].name)
