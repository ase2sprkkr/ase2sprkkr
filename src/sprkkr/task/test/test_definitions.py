if __package__:
   from .init_tests import TestCase
else:
   from init_tests import TestCase, patch_package
   __package__ = patch_package()

import os
from ..tasks import Task

class TestDefinitions(TestCase):

  def test_definitions(self):
      path = os.path.join(os.path.dirname(__file__), '../examples')

      def check(t, name):
          self.assertEqual(t.name, name)
          uname = name.upper()
          self.assertEqual(t.__class__, Task)
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
          self.assertTrue(i[:-3] in Task.definitions())
          td = Task.definitions()[i[:-3]]
          t = td.read_from_file(filename)

          name = i[:-3]
          check(t, name)
          t = Task.from_file(filename)
          check(t, name)
        except Exception as e:
          raise Exception(f'Parsing of "{i}" failed')

