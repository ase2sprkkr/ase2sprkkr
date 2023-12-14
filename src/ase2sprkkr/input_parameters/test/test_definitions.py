import os
import tempfile
import re
import pytest

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..input_parameters import InputParameters


class TestDefinitions(TestCase):

  def change_task(self):
      ip=InputParameters.create_task('DOS')
      ip.ENERGY.EMIN=2
      ip.ENERGY.EMAX=5
      ip.change_task('SCF')
      assert ip.ENERGY.EMIN == 2
      with pytest.raises(AttributeError):
        ip.ENERGY.EMIN

  def jxc(self):
      ip=InputParameters.create_task('SCF')
      ip.MODE.MODE = 'nrel'
      with pytest.warns(UserWarning):
          ip.change_task('jxc')
      with pytest.warns(UserWarning):
          ip.MODE.MODE = 'srel'
      ip.MODE.MODE = 'frel'
      ip.change_task('scf')

  def test_defaults(self):
      for i in InputParameters.definitions:
          ip=InputParameters.create_input_parameters(i)
          df= ip._definition
          ip2 = df.read_from_string(ip.to_string())
          self.assertEqual(ip.to_dict(), ip2.to_dict())
          if i == 'SCF':
              ip.MODE.MDIR[1]=1.,1.,1.
              ip.MODE.MDIR[4]=1.,1.,1.
              ip2 = df.read_from_string(ip.to_string())
              self.assertEqual(ip.to_dict(), ip2.to_dict())

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
          self.assertTrue(i[:-3].upper() in InputParameters.definitions)
          td = InputParameters.definitions[i[:-3].upper()]
          t = td.read_from_file(filename)

          name = i[:-3]
          check(t, name)
          t = InputParameters.from_file(filename)
          check(t, name)

          if name == 'scf':
             self.assertFalse('MODE' in t.to_string())
          if name == 'arpes':
             for i in [ 'AIPES', 'SPLEED', 'BAND' ]:
                 t.TASK.TASK = i

        except Exception as e:
          raise Exception(f'Parsing of "{i}" failed with the reason: \n {e}').with_traceback(e.__traceback__)

      td = InputParameters.definitions['PHAGEN']
      td['TASK']['TASK'].default_value = 'scf'
      t = td.read_from_file(os.path.join(path,'phagen.in'))
      self.assertEqual(t.TASK.TASK(), 'PHAGEN')

      with tempfile.TemporaryFile("w+") as fp:
          t.save_to_file(fp)
          fp.seek(0)
          out = fp.read()
          self.assertTrue(re.sub(r'\s','', out).endswith('TASKPHAGEN'))
      td['TASK']['TASK'].default_value = 'PHAGEN'
