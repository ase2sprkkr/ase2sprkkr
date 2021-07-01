if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
_package__, __name__ = patch_package(__package__, __name__)

import unittest
import pyparsing as pp
from ...common import grammar_types as gt
from .. import task_definitions as cd
from .. import tasks as tasks
from ...common.conf_containers import Section, CustomSection
from ...common.options import Option, CustomOption
import io
import numpy as np
from ...common.grammar import generate_grammar

class TestTask(TestCase):

  def test_section_delimiter_value(self):
     grammar = cd.TaskDefinition.grammar_of_delimiter()
     grammar = 'a' + grammar + 'b'
     for w in ['a b','a\n b','a\n\n b','a\n \n b','a \n\n b', 'a\n\n\n b']:
         self.assertRaises(pp.ParseException, lambda: grammar.parseString(w, True))
     for w in ['a\nb','a \nb','a\n  \nb','a  \n \nb','a \n\t\nb', 'a\n\n\nb']:
         self.assertEqual(['a','b'], grammar.parseString(w, True).asList())


  def test_custom_value(self):
     with generate_grammar():
       cv = cd.SectionDefinition.custom_member_grammar(['aaa'])
     self.assertTrue('\n' not in cv.whiteChars)
     def assertParse(text,result):
         assert cv.parseString(text, True)[0] == result

     def assertNotValid(text):
         self.assertRaises(pp.ParseException, lambda: cv.parseString(text, True))

     assertParse("bbb=1", ('bbb', 1))
     assertParse("bbb=1.3", ('bbb', 1.3))
     assertParse("bbb", ('bbb', True))
     assertNotValid("aaa")
     assertNotValid("aaa=1")
     with generate_grammar():
        cv = cv + 'a'
     assertParse("bbb=1a", ('bbb', 1)) #in this case, 'a' is in result[1]
     assertNotValid("bbb=1\na")



  def test_task_definition(self):
    V = cd.ValueDefinition

    task_def = cd.TaskDefinition.from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        V('NE', gt.SetOf(int, min_length=1)),
        V('Ime', float, 0.0),
      ],
      'SITES' : [
        V('NL', int)
      ]
    })

    def parse(text, _grammar=None):
        _grammar = _grammar or grammar
        return _grammar.parseString(text, True)

    def assertNotValid(text, grammar=None):
      self.assertRaises(pp.ParseBaseException, lambda: parse(text, grammar))

    def assertParse(text, value, _grammar=None):
      out = parse(text, _grammar).asList()
      self.assertEqual(len(out),1)
      self.assertEqual(out[0], value)

    def ar(x):
      return np.atleast_1d(x)

    grammar = task_def.sections['SITES'].values['NL'].grammar()
    assertParse("NL=3", ('NL', 3))
    grammar = task_def.sections['ENERGY'].values['NE'].grammar()
    assertParse("NE={3}", ('NE', ar(3)))
    grammar = task_def.sections['ENERGY'].values['Ime'].grammar()
    assertParse("Ime= 0.5", ('Ime', 0.5))
    grammar = task_def.sections['ENERGY'].values['GRID'].grammar()
    assertParse("GRID={3}", ('GRID', ar(3)))
    grammar = task_def.sections['ENERGY'].grammar()
    assertParse("ENERGY Ime= 0.5", ('ENERGY', {'Ime':0.5}))
    assertParse("ENERGY Ime= 0.5 NE={5}",('ENERGY', {'Ime':0.5, 'NE':ar(5)}) )
    assertParse("""ENERGY Ime= 0.5
                                   NE={5}""", ('ENERGY', {'Ime':0.5, 'NE':ar(5)}) )
    assertParse("""ENERGY Ime= 0.5

                                   NE={5}""", ('ENERGY', {'Ime':0.5, 'NE':ar(5)}) )
    assertNotValid("""ENERGY Ime= 0.5

NE={5}""")

    #fixed value = 3
    assertNotValid("""ENERGY GRID={1}""")
    #no space before a section name
    assertNotValid(""" ENERGY GRID={3}""")
    assertParse("""ENERGY GRID={3}""", ('ENERGY', {'GRID':ar(3)}))

    grammar = task_def.grammar()
    assertParse("""ENERGY GRID={3}""", {'ENERGY': {'GRID':ar(3)}})
    assertParse("""ENERGY GRID={3}
                     NE={300}


          """, {'ENERGY': {'GRID':3, 'NE':300}})
    assertParse("""ENERGY
                     NE={300}
                     GRID={3}


SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300}, 'SITES':{'NL':2}} )

    #custom values
    with generate_grammar():
      grammar = task_def['ENERGY']._values_grammar()
    assertParse("""GRID={3}
                     NE={300}
                     """, {'GRID':ar(3), 'NE':ar(300)})

    assertParse("""GRID={3}
                   NE={300}
                   NXXX=5""", {'GRID':ar(3), 'NE':ar(300), 'NXXX': 5})

    grammar = task_def.grammar()
    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX=5


SITES NL=2""", {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': 5},
                              'SITES':{'NL':2}}
    )


    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX

SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True},
                              'SITES':{'NL':2}}
    )

    #SITES do not start on the begin of the line, so it is not the start of the section
    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX

     SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True,
                              'SITES': True, 'NL':2 }}
    )


    assertParse(""" ENERGY GRID={3}
                     NE={300}
                     NXXX

  SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True, 'SITES': True, 'NL': 2}}
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

      {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3}
      })

    #multiline custom section
    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

XSITES NR=3 NF=1
                     NZ=5.5
              """,

      {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3, 'NF':1,'NZ':5.5}
      })

    #custom section with a flag
    assertParse("""ENERGY GRID={3}
                     NXXX
                     NZZZ=4
                     NE={300}


SITES NL=2
 
XSITES NR=3 FLAG
                     FLOAT=3.5
              """,

         {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True, 'NZZZ':4},
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
    self.assertEqual(task.find('NL').get_path(), 'SITES.NL')
    task.find('NL').set(3)
    self.assertEqual(task['SITES'].NL(), 3)
    self.assertTrue(isinstance(task['XSITES'], CustomSection))

    output = io.StringIO()
    task.save_to_file(output)
    output.seek(0)
    task2 = task_def.read_from_file(output)
    self.assertEqual(str(task.to_dict()), str(task2.to_dict()))
