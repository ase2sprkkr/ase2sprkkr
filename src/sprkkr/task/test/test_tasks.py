from init_tests import patch_package
__package__ = patch_package()

print(__package__)
import unittest
import pyparsing
from ...common import grammar_types as gt
from .. import task_definitions as cd

class OptionsTest(unittest.TestCase):

  def test_custom_value(self):
     cv = cd.custom_value(['aaa'])

     def assertParse(text,result):
         assert cv.parseString(text, True)[0] == result

     def assertNotValid(text):
         self.assertRaises(pp.ParseException, lambda: cv.parseString(text, True))
         
     assertParse("bbb=1.3", ('bbb', 1.3))
     assertParse("bbb=1", ('bbb', 1))
     assertParse("bbb", ('bbb', True))
     assertNotValid("aaa")
     assertNotValid("aaa=1")
     import sys
     sys.exit()


  def test_task_definition(self):
    V = cd.ValueDefinition

    conf_def = cd.TaskDefinition.from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        V('NE', gt.SetOf(int, length=1)),
        V('Ime', float, 0.0),
      ],
      'SITES' : [
        V('NL', int)
      ]
    })

    def parse(text):
        return grammar.parseString(text, True)

    def assertNotValid(text):
      self.assertRaises(pyparsing.ParseException, lambda: parse(text))

    def assertParse(text, value):
      out = parse(text).asList()
      self.assertEqual(len(out),1)
      self.assertEqual(out[0], value)

    grammar = conf_def.sections['SITES'].values['NL'].grammar()
    assertParse("NL=3", ('NL', 3))
    grammar = conf_def.sections['ENERGY'].values['NE'].grammar()
    assertParse("NE={3}", ('NE', 3))
    grammar = conf_def.sections['ENERGY'].values['Ime'].grammar()
    assertParse("Ime= 0.5", ('Ime', 0.5))
    grammar = conf_def.sections['ENERGY'].values['GRID'].grammar()
    assertParse("GRID={3}", ('GRID', 3))
    grammar = conf_def.sections['ENERGY'].grammar()
    assertParse("ENERGY Ime= 0.5", ('ENERGY', {'Ime':0.5}))
    assertParse("ENERGY Ime= 0.5 NE={5}",('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertParse("""ENERGY Ime= 0.5
                                   NE={5}""", ('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertNotValid("""ENERGY Ime= 0.5

                                   NE={5}""")
    assertNotValid(""" ENERGY GRID={1}""")
    assertParse(""" ENERGY GRID={3}""", ('ENERGY', {'GRID':3}))

    grammar = conf_def.grammar()
    assertParse(""" ENERGY GRID={3}""", {'ENERGY': {'GRID':3}})
    assertParse(""" ENERGY GRID={3}
                     NE={300}


          """, {'ENERGY': {'GRID':3, 'NE':300}})
    assertParse(""" ENERGY GRID={3}
                     NE={300}


              SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300}, 'SITES':{'NL':2}} )

    #custom values

    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300}, 'SITES':{'NL':2}} )
