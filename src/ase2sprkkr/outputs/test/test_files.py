if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from os import listdir
from os.path import dirname, join as path_join, isfile
import os.path
from ..output_files import OutputFile

class TestOutput(TestCase):

  def test_output(self):
    dire = path_join(dirname(dirname(__file__)), 'output_files', 'examples')
    for i in listdir(dire):
        fname = path_join(dire, i)
        if isfile(fname):
            out=OutputFile.from_file(fname)
            ext=fname.rsplit('.',1)[1]
            if ext=='spc':
               self.assertEqual('ARPES', out.KEYWORD())
               self.assertEqual((200,160), out.ENERGY().shape)
