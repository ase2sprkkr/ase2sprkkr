if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


import os, glob, re  # NOQA


class TestCommon(TestCase):

  def test(self):
      dire = os.path.join( os.path.join(os.path.dirname(__file__), '..') )
      examples = glob.glob( os.path.join(dire, 'examples', '*', 'run.sh') )
      cwd = os.getcwd()
      try:
        for i in examples:
          with open(i, 'r') as f:
               os.chdir(os.path.dirname(i))
               lines = f.readlines()
               cmd = lines[1]
               cmd = re.sub(' -a( |$)','\\1', cmd)
               self.assertEqual(0, os.system(cmd))
      finally:
        os.chdir(cwd)
