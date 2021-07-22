if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..input_parameters import InputParameters

class TestDefinitions(TestCase):

  def test_definitions(self):
      path = os.path.join(os.path.dirname(__file__), '../examples')

      def check(t, name):
          self.assertEqual(t.name, name)
          uname = name.upper()
          self.assertEqual(t.__class__, InputParameters)
          self.assertEqual(t['CONTROL']['ADSI'](), uname)
          self.assertEqual(name != 'scf', uname in t['TASK'])
          if name != 'scf':
            self.assertEqual(t['TASK'][uname](), True)

      #for i in os.listdir( path ):
      for i in ['arpes.in']:
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

