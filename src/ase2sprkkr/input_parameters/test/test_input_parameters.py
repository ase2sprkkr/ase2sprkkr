if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import unittest
import pyparsing as pp
from ...common import grammar_types as gt
from .. import input_parameters_definitions as cd
from .. import input_parameters as input_parameters
from ...common.configuration_containers import Section, CustomSection
from ...common.options import Option, CustomOption
import io, re
import numpy as np
from ...common.grammar import generate_grammar
import pytest

class TestInputParameters(TestCase):

  def test_section_delimiter_value(self):
     grammar = cd.InputParametersDefinition.grammar_of_delimiter()
     grammar = 'a' + grammar + 'b'
     for w in ['a b','a\n b','a\n\n b','a\n \n b','a \n\n b', 'a\n\n\n b']:
         self.assertRaises(pp.ParseException, lambda: grammar.parseString(w, True))
     for w in ['a\nb','a \nb','a\n  \nb','a  \n \nb','a \n\t\nb', 'a\n\n\nb']:
         self.assertEqual(['a','b'], grammar.parseString(w, True).asList())


  def test_custom_value(self):
     with generate_grammar():
       cv = cd.InputSectionDefinition.custom_member_grammar(['aaa'])
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

  def test_dangerous_value(self):
     with generate_grammar():
       cv = cd.InputValueDefinition('aaa', int, required=True).grammar()
     def assertParse(text,result):
         assert cv.parseString(text, True)[0] == result
     def assertParseDangerous(text,result):
         out = cv.parseString(text, True)[0]
         assert out[0] == result[0]
         self.assertEqual(out[1].value, result[1])
     def assertNotValid(text):
         self.assertRaises(pp.ParseBaseException, lambda: cv.parseString(text, True))

     assertParse("aaa=1", ('aaa', 1))
     assertNotValid("aaa")
     assertNotValid("aaa=sss")
     assertNotValid("aaa=1.0")

     with generate_grammar():
       cv = cd.InputValueDefinition('aaa', int, required=True).grammar(allow_dangerous=True)

     assertParse("aaa=1", ('aaa', 1))
     assertParseDangerous("aaa=1.0", ('aaa', 1.0))
     assertParseDangerous("aaa=sss", ('aaa', 'sss'))
     assertParseDangerous("aaa={2,3}", ('aaa', np.array([2,3])))

  def test_write_condition(self):
    V = cd.InputValueDefinition

    input_parameters_def = cd.InputParametersDefinition.from_dict({
      'ENERGY' : [
        V('E', 1)
      ]})
    del input_parameters_def._members['TASK']
    id=input_parameters_def.read_from_file(io.StringIO('ENERGY E=1'))
    self.assertEqual(id.to_string(), "ENERGY\n\tE=1\n")
    input_parameters_def['ENERGY']['E'].write_condition = lambda o: True
    self.assertEqual(id.to_string(), "ENERGY\n\tE=1\n")
    input_parameters_def['ENERGY']['E'].write_condition = lambda o: False
    self.assertEqual(id.to_string(), "ENERGY\n")
    input_parameters_def['ENERGY'].write_condition = lambda o: True
    self.assertEqual(id.to_string(), "ENERGY\n")
    input_parameters_def['ENERGY'].write_condition = lambda o: False
    self.assertEqual(id.to_string(), "")

  def test_input_parameters_definition(self):
    V = cd.InputValueDefinition

    input_parameters_def = cd.InputParametersDefinition.from_dict({
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

    grammar = input_parameters_def.sections['SITES'].values['NL'].grammar()
    assertParse("NL=3", ('NL', 3))
    grammar = input_parameters_def.sections['ENERGY'].values['NE'].grammar()
    assertParse("NE={3}", ('NE', ar(3)))
    grammar = input_parameters_def.sections['ENERGY'].values['Ime'].grammar()
    assertParse("Ime= 0.5", ('Ime', 0.5))
    grammar = input_parameters_def.sections['ENERGY'].values['GRID'].grammar()
    assertParse("GRID={3}", ('GRID', ar(3)))
    grammar = input_parameters_def.sections['ENERGY'].grammar()
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

    grammar = input_parameters_def.grammar()
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
      grammar = input_parameters_def['ENERGY']._values_grammar()
    assertParse("""GRID={3}
                     NE={300}
                     """, {'GRID':ar(3), 'NE':ar(300)})

    assertParse("""GRID={3}
                   NE={300}
                   NXXX=5""", {'GRID':ar(3), 'NE':ar(300), 'NXXX': 5})

    grammar = input_parameters_def.grammar()
    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX=5


SITES NL=2""", {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': 5},
                              'SITES':{'NL':2}}
    )


    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX

SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True},
                              'SITES':{'NL':2}}
    )

    #SITES do not start on the begin of the line, so it is not the start of the section
    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX

     SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True,
                              'SITES': True, 'NL':2 }}
    )


    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX

  SITES NL=2

              """, {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True, 'SITES': True, 'NL': 2}}
    )




    assertNotValid("""ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

SITES NL=3
              """)

    #custom section
    assertParse("""ENERGY GRID={3}
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

    ips=input_parameters_def.read_from_file(io.StringIO("""ENERGY GRID={3}
                     NE={300,200}
                     NXXX


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """))
    self.assertTrue(isinstance(ips, input_parameters.InputParameters))
    self.assertTrue(isinstance(ips['ENERGY'], Section))
    self.assertTrue(isinstance(ips['ENERGY']['NE'], Option))
    self.assertTrue(isinstance(ips['ENERGY']['NXXX'], CustomOption))
    self.assertTrue(isinstance(ips.ENERGY.NXXX, CustomOption))
    ips.ENERGY.NXXX.set(5)
    self.assertEqual(ips.ENERGY.NXXX(), 5)
    ips.ENERGY.NXXX.remove()
    self.assertFalse(hasattr(ips.ENERGY, 'NXXX'))

    self.assertEqual(ips['ENERGY'].NE(), np.array((300,200)))
    cps = ips.copy(copy_values=True)
    cps.ENERGY.NE[0]=100
    self.assertEqual(ips['ENERGY'].NE(), np.array((300,200)))
    dps = cps.copy()
    dps.ENERGY.NE[0]=150
    self.assertEqual(cps['ENERGY'].NE(), np.array((150,200)))
    cps.ENERGY.NE[0]=180
    self.assertEqual(dps['ENERGY'].NE(), np.array((180,200)))

    self.assertEqual(ips['SITES'].NL(), 2)
    self.assertEqual(ips.find('NL').get_path(), 'SITES.NL')
    ips.find('NL').set(3)
    self.assertEqual(ips['SITES'].NL(), 3)
    self.assertTrue(isinstance(ips['XSITES'], CustomSection))

    output = io.StringIO()
    ips.save_to_file(output)
    output.seek(0)
    ips2 = input_parameters_def.read_from_file(output)
    self.assertEqual(str(ips.as_dict()), str(ips2.as_dict()))

    with pytest.raises(ValueError):
        ips.ENERGY.NE = 'sss'
    ips.ENERGY.NE.set_dangerous('sss')
    output = io.StringIO()
    ips.save_to_file(output)
    output.seek(0)
    with pytest.raises(pp.ParseBaseException):
        ips2 = input_parameters_def.read_from_file(output)
    output.seek(0)
    ips2 = input_parameters_def.read_from_file(output, allow_dangerous=True)
    self.assertEqual('sss', ips.ENERGY.NE())
    self.assertEqual(str(ips.as_dict()), str(ips2.as_dict()))

    #numbered_arrays
    input_parameters_def = cd.InputParametersDefinition.from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        V('NE', gt.SetOf(int, min_length=1)),
        V('Ime', float, 0.0),
      ],
      'SITES' : [
        V('NL', int)
      ]
    })

    assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0");
    assertParse("ENERGY NE=1 Ime=0.5", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : 0.5}});
    input_parameters_def.sections['ENERGY']['Ime'].is_numbered_array = True
    with generate_grammar():
      grammar = input_parameters_def._values_grammar()

    assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0");
    assertParse("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : {'def': 0.5,  1:0.4, 5:0.8} }, });

    ip=input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
    self.assertEqual(0.5, ip.ENERGY.Ime())
    self.assertEqual(0.8, ip.ENERGY.Ime[5])
    ip.ENERGY.Ime[1] = 0.7
    ip.ENERGY.Ime[9] = 0.9
    self.assertEqual('ENERGY GRID={3} NE={1} Ime=0.5 Ime1=0.7 Ime5=0.8 Ime9=0.9', re.sub(r'\s+',' ', ip.ENERGY.to_string()).strip() )
    ip.ENERGY.Ime[[1,4,9]] = 0.2
    self.assertEqual([0.2,0.8,0.2], ip.ENERGY.Ime[1,5,9])
    ip.ENERGY.Ime[[1,4,9]] = 0.3
    self.assertEqual([0.3,0.8,0.3], ip.ENERGY.Ime[[1,5,9]])
    ip.ENERGY.Ime[5:9] = 0.2
    self.assertEqual([0.2,0.2,0.2,0.2,0.3], ip.ENERGY.Ime[5:10])
    self.assertRaises(KeyError, lambda: ip.ENERGY.Ime[5.0])
    def e():
       ip.ENERGY.Ime[5.0]=1
    self.assertRaises(KeyError, e)
    self.assertRaises(KeyError, lambda: ip.ENERGY.Ime['5'])
    def e():
       ip.ENERGY.Ime['5']=1
    self.assertRaises(KeyError, e)

    ip=input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
    with pytest.raises(ValueError):
      ip.ENERGY.Ime = 'ss'
    ip.ENERGY.Ime.set_dangerous('uu')
    self.assertEqual(ip.ENERGY.Ime(),'uu')
    out=ip.ENERGY.to_string()
    self.assertEqual('ENERGY GRID={3} NE={1} Ime=uu Ime1=0.4 Ime5=0.8', re.sub(r'\s+',' ', out).strip() )
    with pytest.raises(pp.ParseBaseException):
         input_parameters_def.read_from_file(io.StringIO(out))
    ip=input_parameters_def.read_from_file(io.StringIO(out), allow_dangerous=True)
    self.assertEqual(ip.ENERGY.Ime(all_values=True), {'def': 'uu', 1:0.4, 5:0.8 } )
    ip.ENERGY.Ime = 1.0

    with pytest.raises(ValueError):
      ip.ENERGY.Ime[5] = 'ss'
    ip.ENERGY.Ime.set_dangerous('yy', index=5)
    out=ip.ENERGY.to_string()
    self.assertEqual('ENERGY GRID={3} NE={1} Ime=1.0 Ime1=0.4 Ime5=yy', re.sub(r'\s+',' ', ip.ENERGY.to_string()).strip() )
    with pytest.raises(pp.ParseBaseException):
         input_parameters_def.read_from_file(io.StringIO(out))
    ip=input_parameters_def.read_from_file(io.StringIO(out), allow_dangerous=True)
    self.assertEqual(ip.ENERGY.Ime(all_values=True), {'def': 1.0, 1:0.4, 5:'yy' } )
