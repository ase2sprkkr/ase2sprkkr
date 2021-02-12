from init_tests import patch_package
__package__ = patch_package()
import unittest
import pyparsing
from .. import sprkkr_io_types as sit
from .. import sprkkr_io_entry_def as sied

class IosTest(unittest.TestCase):

  def test_types(self):
    Error = "ERROR"

    type = sit.Integer()

    def test(val, res):
        try:
          out = type.grammar().parseString(val, True)
          out = type.read(out)
          assert out == res, "{} should be {} and is {} for type {}".format(val, res, out, type.__class__.__name__)
        except (ValueError, pyparsing.ParseException):
          assert res is Error, "{} should be validated as {} for type {}".format(val, res, type.__class__.__name__)
          return
        assert type.validate(out)

    for val, res in [
        ('1', 1),
        ('a', Error),
        ('1.1', Error),
        ('1.1 a', Error),
        ('-1', -1)
        ]:
        test(val, res)

    type = sit.Real()
    for val, res in [
        ('1', 1.0),
        ('a', Error),
        ('1.1', 1.1),
        ('1.1 a', Error),
        ('-1e-2', -1e-2)
        ]:
        test(val, res) 

    type = sit.String()
    for val, res in [
         ('a', 'a'),
         ('a a', Error),
         ('1','1')
         ]:
         test(val, res)

    type = sit.SetOf(int)
    for val, res in [
        ('a{}', Error),
        ('{a}', Error),
        ('{1,2}', [1,2]),
        ('{1,a}', Error)
        ]:
        test(val, res)

  def test_entry_def(self):
    V = sied.EntryValueDef

    entry_def = sied.EntryDef.from_dict({
      'ENERGY' : [
        V('GRID', sit.SetOf(int, length=1), fixed_value=3),
        V('NE', sit.SetOf(int, length=1)),
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

    grammar = entry_def.sections['SITES'].values['NL'].grammar()
    assertParse("NL=3", ('NL', 3))
    grammar = entry_def.sections['ENERGY'].values['NE'].grammar()
    assertParse("NE={3}", ('NE', 3))
    grammar = entry_def.sections['ENERGY'].values['Ime'].grammar()
    assertParse("Ime= 0.5", ('Ime', 0.5))
    grammar = entry_def.sections['ENERGY'].values['GRID'].grammar()
    assertParse("GRID={3}", ('GRID', 3))
    grammar = entry_def.sections['ENERGY'].grammar()
    assertParse("ENERGY Ime= 0.5", ('ENERGY', {'Ime':0.5}))
    assertParse("ENERGY Ime= 0.5 NE={5}",('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertParse("""ENERGY Ime= 0.5
                                   NE={5}""", ('ENERGY', {'Ime':0.5, 'NE':5}) )
    assertNotValid("""ENERGY Ime= 0.5

                                   NE={5}""")
    assertNotValid(""" ENERGY GRID={1}""")
    assertParse(""" ENERGY GRID={3}""", ('ENERGY', {'GRID':3}))

    grammar = entry_def.grammar()
    assertParse(""" ENERGY GRID={3}""", {'ENERGY': {'GRID':3}})
    assertParse(""" ENERGY GRID={3}
                     NE={300}


          """, {'ENERGY': {'GRID':3, 'NE':300}})
    assertParse(""" ENERGY GRID={3}
                     NE={300}


              SITES NL=2""", {'ENERGY': {'GRID':3, 'NE':300}, 'SITES':{'NL':2}} )








    

       
       
     
  
