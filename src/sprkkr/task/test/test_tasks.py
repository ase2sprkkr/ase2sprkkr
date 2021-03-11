if __package__:
   from .init_tests import TestCase
else:
   from init_tests import TestCase, patch_package
   __package__ = patch_package()

import unittest
import pyparsing as pp
from ...common import grammar_types as gt
from .. import task_definitions as cd
from .. import tasks as tasks
from ...common.conf_containers import Section, CustomSection
from ...common.options import Option, CustomOption
import io
import numpy as np

class TestTask(TestCase):

  def test_custom_value(self):

     cv = cd.SectionDefinition.custom_value(['aaa'])
     def assertParse(text,result):
         assert cv.parseString(text, True)[0] == result

     def assertNotValid(text):
         self.assertRaises(pp.ParseException, lambda: cv.parseString(text, True))

     assertParse("bbb=1", ('bbb', 1))
     assertParse("bbb=1.3", ('bbb', 1.3))
     assertParse("bbb", ('bbb', True))
     assertNotValid("aaa")
     assertNotValid("aaa=1")


  def test_task_definition(self):
    V = cd.ValueDefinition

    task_def = cd.TaskDefinition.from_dict({
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
      self.assertRaises(pp.ParseException, lambda: parse(text))

    def assertParse(text, value):
      out = parse(text).asList()
      self.assertEqual(len(out),1)
      self.assertEqual(out[0], value)

    grammar = task_def.sections['SITES'].values['NL'].grammar()
    assertParse("NL=3", ('NL', 3))
    grammar = task_def.sections['ENERGY'].values['NE'].grammar()
    assertParse("NE={3}", ('NE', 3))
    grammar = task_def.sections['ENERGY'].values['Ime'].grammar()
    assertParse("Ime= 0.5", ('Ime', 0.5))
    grammar = task_def.sections['ENERGY'].values['GRID'].grammar()
    assertParse("GRID={3}", ('GRID', 3))
    grammar = task_def.sections['ENERGY'].grammar()
    assertParse("ENERGY Ime= 0.5", ('ENERGY', {'Ime':0.5}))
    assertParse("ENERGY Ime= 0.5 NE={5}",('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertParse("""ENERGY Ime= 0.5
                                   NE={5}""", ('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertNotValid("""ENERGY Ime= 0.5

                                   NE={5}""")
    assertNotValid(""" ENERGY GRID={1}""")
    assertParse(""" ENERGY GRID={3}""", ('ENERGY', {'GRID':3}))

    grammar = task_def.grammar()
    assertParse(""" ENERGY GRID={3}""", {'ENERGY': {'GRID':3}})
    assertParse(""" ENERGY GRID={3}
                     NE={300}


          """, {'ENERGY': {'GRID':3, 'NE':300}})
    assertParse(""" ENERGY
                     GRID={3}
                     NE={300}


              SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300}, 'SITES':{'NL':2}} )

    #custom values and sections

    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300, 'NXXX': True},
                              'SITES':{'NL':2}}
    )


    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX

              SITES NL=2

              """, {'ENERGY': {'GRID':3, 'NE':300, 'NXXX': True},
                              'SITES':{'NL':2}}
    )



    assertNotValid(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2

              SITES NL=3
              """)

    #custom section
    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2

              XSITES NR=3
              """,

      {'ENERGY': {'GRID':3, 'NE':300, 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3}
      })

    #multiline custom section
    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2

              XSITES NR=3 NF=1
                     NZ=5.5
              """,

      {'ENERGY': {'GRID':3, 'NE':300, 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3, 'NF':1,'NZ':5.5}
      })

    #custom section with a flag
    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX


              SITES NL=2

              XSITES NR=3 FLAG
                     FLOAT=3.5
              """,

      {'ENERGY': {'GRID':3, 'NE':300, 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3, 'FLAG' : True, 'FLOAT': 3.5}
      })

    task=task_def.read_from_file(io.StringIO(""" ENERGY GRID={3}
                     NE={300,200}
                     NXXX


              SITES NL=2

              XSITES NR=3 FLAG
                     FLOAT=3.5
              """))
    self.assertTrue(isinstance(task, tasks.Task))
    self.assertTrue(isinstance(task['ENERGY'], Section))
    self.assertTrue(isinstance(task['ENERGY']['NE'], Option))
    self.assertTrue(isinstance(task['ENERGY']['NXXX'], CustomOption))
    self.assertEqual(task['ENERGY'].NE(), np.array((300,200)))
    self.assertEqual(task['SITES'].NL(), 2)
    self.assertTrue(isinstance(task['XSITES'], CustomSection))

    output = io.StringIO()
    task.save_to_file(output)
    output.seek(0)
    task2 = task_def.read_from_file(output)
    self.assertEqual(str(task.to_dict()), str(task2.to_dict()))
