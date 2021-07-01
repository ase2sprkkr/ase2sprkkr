if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..misc import add_to_signature
import inspect
import tempfile
import asyncio
from ..process_output_reader import AsyncioFileReader

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

  def test_asyncio_file_reader(self):
    with tempfile.TemporaryFile(mode='w+b') as tfile:
      tfile.write(b"""First line\nSecond line\nHi, this is a very long file!! really!!! 12345123456789""")
      tfile.seek(0)
      ar = AsyncioFileReader(tfile, buffersize=5)
      self.assertAsyncEqual(b"First line\n", ar.readline())
      self.assertAsyncEqual(b"Second line\n", ar.readline())
      self.assertAsyncEqual(b"Hi, this", ar.readuntil(b'is'))
      self.assertAsyncEqual(b" is", ar.readuntil(b'is'))
      self.assertAsyncEqual(b" a very long file!! really!!!", ar.readuntil(b'!!!'))
      self.assertAsyncEqual(b" 12345123456", ar.readuntil(b'123456'))
      with self.assertRaises(asyncio.IncompleteReadError):
        self.runAsync(ar.readuntil(b'123456'))
