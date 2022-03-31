if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..input_parameters import InputParameters
import tempfile
import re

class TestDefinitions(TestCase):

  def test_definitions(self):
      path = os.path.join(os.path.dirname(__file__), '../examples')

      def check(t, name):
          self.assertEqual(t.name, name)
          uname = name.upper()
          self.assertEqual(t.__class__, InputParameters)
          self.assertEqual(t['CONTROL']['ADSI'](), uname)
          self.assertEqual(t.TASK.TASK(), uname)

      for i in os.listdir( path ):
        try:
          if not i.endswith('.in'):
             continue
          filename = os.path.join(path, i)
          self.assertTrue(i[:-3].upper() in InputParameters.definitions())
          td = InputParameters.definitions()[i[:-3].upper()]
          t = td.read_from_file(filename)

          name = i[:-3]
          check(t, name)
          t = InputParameters.from_file(filename)
          check(t, name)
        except Exception as e:
          raise Exception(f'Parsing of "{i}" failed with the reason: \n {e}').with_traceback(e.__traceback__)

      td = InputParameters.definitions()['PHAGEN']
      td['TASK']['TASK'].default_value = 'scf'
      t = td.read_from_file(os.path.join(path,'phagen.in'))
      self.assertEqual(t.TASK.TASK(), 'PHAGEN')

      with tempfile.TemporaryFile("w+") as fp:
          t.save_to_file(fp)
          fp.seek(0)
          out = fp.read()
          self.assertTrue(re.sub(r'\s','', out).endswith('TASKPHAGEN'))
      td['TASK']['TASK'].default_value = 'PHAGEN'
