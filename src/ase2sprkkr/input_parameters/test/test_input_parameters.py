import pyparsing as pp
import io
import re
import numpy as np
import pytest
from functools import partial

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:  # Just a linter worshiping
    from ...common.grammar import generate_grammar
    from ...common import grammar_types as gt
    from .. import input_parameters_definitions as cd
    from .. import input_parameters as input_parameters
    from ...common.configuration_containers import Section, CustomSection
    from ...common.options import Option, CustomOption
    from ...common.configuration_definitions import gather, switch
    from ...common.generated_configuration_definitions import Length
    from ...common.warnings import DataValidityError

V = cd.InputValueDefinition


def ar(x):
   return np.atleast_1d(x)


class TestInputParameters(TestCase):

  def assertParse(self, text, value, grammar):
      out = self.parse(text, grammar).asList()
      self.assertEqual(len(out),1)
      out = out[0]
      if hasattr(out, 'to_dict'):
          out = out.to_dict()
      self.assertEqual(out, value)

  def parse(self, text, grammar):
        if not isinstance(grammar, pp.ParserElement):
           grammar = grammar()
        return grammar.parseString(text, True)

  def assertNotValid(self, text, grammar=None):
      self.assertRaises(pp.ParseBaseException, lambda: self.parse(text, grammar))

  def test_section_delimiter_value(self):
     grammar = cd.InputParametersDefinition.grammar_of_delimiter()
     grammar = 'a' + grammar + 'b'
     for w in ['a b','a\n b','a\n\n b','a\n \n b','a \n\n b', 'a\n\n\n b']:
         self.assertRaises(pp.ParseException, lambda: grammar.parseString(w, True))
     for w in ['a\nb','a \nb','a\n  \nb','a  \n \nb','a \n\t\nb', 'a\n\n\nb']:
         self.assertEqual(['a','b'], grammar.parseString(w, True).asList())
#

  def test_custom_value(self):
     with generate_grammar():
       cv = cd.InputSectionDefinition.custom_member_grammar(lambda x: x[0] not in ['aaa'])
     self.assertTrue('\n' not in cv.whiteChars)

     assertParse = partial(self.assertParse, grammar=lambda: cv)
     assertNotValid = partial(self.assertNotValid, grammar=lambda: cv)

     assertParse("bbb=1", ('bbb', 1))
     assertParse("bbb=1.3", ('bbb', 1.3))
     assertParse("bbb", ('bbb', True))
     assertNotValid("aaa")
     assertNotValid("aaa=1")
     with generate_grammar():
        cv = cv + 'a'
     out = self.parse("bbb=1 a", cv)
     self.assertEqual(list(out), [('bbb', 1), 'a'])
     out = self.parse("bbb=1a a", cv)
     self.assertEqual(list(out), [('bbb', "1a"), 'a'])
     assertNotValid("bbb=1\na")

  def test_dangerous_value(self):
     with generate_grammar():
       cv = cd.InputValueDefinition('aaa', int, is_required=True).grammar()

     assertParse = partial(self.assertParse, grammar=lambda: cv)
     assertNotValid = partial(self.assertNotValid, grammar=lambda: cv)

     def assertParseDangerous(text,result):
         out = cv.parseString(text, True)[0]
         assert out[0] == result[0]
         self.assertEqual(out[1].value, result[1])

     assertParse("aaa=1", ('aaa', 1))
     assertNotValid("aaa")
     assertNotValid("aaa=sss")
     assertNotValid("aaa=1.0")

     with generate_grammar():
       cv = cd.InputValueDefinition('aaa', int, is_required=True).grammar(allow_dangerous=True)

     assertParse("aaa=1", ('aaa', 1))
     assertParseDangerous("aaa=1.0", ('aaa', 1.0))
     assertParseDangerous("aaa=sss", ('aaa', 'sss'))
     assertParseDangerous("aaa={2,3}", ('aaa', np.array([2,3])))

     vd=cd.InputValueDefinition('aaa', 2, is_required=True)
     o=vd.create_object()
     o.set_dangerous("AAA")
     self.assertEqual(o.to_string(), '\taaa=AAA')
     o.set_dangerous(None)
     self.assertEqual(o.to_string(), '')
     #

  def test_is_required(self):
      ipd = cd.InputParametersDefinition.definition_from_dict({
        'ENERGY' : [
          V('E', 1),
          V('F', int, None, is_required = lambda o: o._container.E()!=2)
        ]
      })

      ip = ipd.create_object()
      with pytest.raises(DataValidityError):
          ip.validate()
      e = ip.ENERGY
      e.F = 2
      ip.validate()
      with pytest.raises(DataValidityError):
          e.F = None
      e.E = 2
      ip.validate()
      e.F = None
      ip.validate()

      assert ipd.read_from_string("ENERGY E=2").ENERGY.to_dict() == { 'E': 2 }
      assert ipd.read_from_string("ENERGY E=1 F=2").ENERGY.to_dict() == { 'E': 1, 'F': 2 }
      with pytest.warns(DataValidityError):
          ipd.read_from_string("ENERGY E=1")

  def test_write_condition(self):
    input_parameters_def = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('E', 1)
      ]})
    # del input_parameters_def._members['TASK']
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
    #

  def test_input_parameters_definition(self):
    input_parameters_def = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        V('NE', gt.SetOf(int, min_length=1)),
        V('Ime', float, 0.0),
      ],
      'SITES' : [
        V('NL', int)
      ]
    })
    grammar = input_parameters_def.sections['SITES'].values['NL'].grammar()
    assertParse = partial(self.assertParse, grammar=lambda: grammar)
    assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

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

    # fixed value = 3
    assertNotValid("""ENERGY GRID={1}""")
    # no space before a section name
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


SITES NL=2""", {'ENERGY': {'NE':300, 'GRID':3}, 'SITES':{'NL':2}} )

    # custom values
    with generate_grammar():
      grammar = input_parameters_def['ENERGY']._grammar_of_values()
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

    # SITES do not start on the begin of the line, so it is not the start of the section
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

    # custom section
    assertParse("""ENERGY GRID={3}
                     NE={300}
                     NXXX

    #numbered_arrays

SITES NL=2

XSITES NR=3
              """,

      {'ENERGY': {'GRID':ar(3), 'NE':ar(300), 'NXXX': True},
                              'SITES':{'NL':2},
                              'XSITES':{'NR':3}
      })

    # multiline custom section
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

    # a custom section with a flag
    assertParse("""ENERGY GRID={3}
                     NXXX
                     NE={300}
                     NZZZ=4


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """,

     {'ENERGY': {'GRID':ar(3), 'NXXX':True, 'NE':ar(300), 'NZZZ':4},
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

    with pytest.raises(DataValidityError):
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
#

  def test_set_values(self):
    input_parameters_def = cd.InputParametersDefinition.definition_from_dict({
        'ENERGY' : [
          V('GRID', gt.SetOf(int, length=1), fixed_value=3),
          V('NE', gt.SetOf(int, min_length=1)),
          V('Ime', float, 0.0),
          V('ENERGY', float, 0.0),
        ],
        'SITES' : [
          V('NL', 1)
        ]
      })
    ips=input_parameters_def.create_object()

    ips.set({'ENERGY.Ime' : 5., 'NE' : 10, 'ENERGY': 7. })
    self.assertEqual(ips.ENERGY.Ime(), 5.)
    self.assertEqual(ips.ENERGY.NE(), 10)
    self.assertEqual(ips.ENERGY.ENERGY(), 7.)
    ips.set({'ENERGY.ENERGY': 4. })
    self.assertEqual(ips.ENERGY.ENERGY(), 4.)
    ips.set('ENERGY.ENERGY', 5. )
    self.assertEqual(ips.ENERGY.ENERGY(), 5.)
    ips.set('ENERGY', 6. )
    self.assertEqual(ips.ENERGY.ENERGY(), 6.)
    ips.set('ENERGY.ENERGYY', 7., unknown='ignore')
    self.assertFalse('ENERGYY' in ips.ENERGY)
    self.assertRaises(KeyError, lambda: ips.set('ENERGY.ENERGYY', 7., unknown='fail'))
    ips.set('ENERGY.ENERGYY', 7., unknown='add')
    self.assertEqual(ips.ENERGY.ENERGYY(), 7.)
#

  def test_repeated_value(self):

    t=lambda x:re.sub(r'\s+',' ',x).strip()
    assertParse = self.assertParse
    assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('B', 2, is_repeated = True),
        V('C', 3) ]})

    data="""ENERGY A=3
    B=2
    C=77"""
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [2], 'C':77}}, ipd.grammar())
    ip=ipd.read_from_string(data)
    self.assertEqual({'ENERGY': {'A' : 3, 'B' : [2], 'C':77}}, ip.to_dict())
    self.assertEqual('ENERGY A=3 B=2 C=77', t(ip.ENERGY.to_string()) )

    data="""ENERGY A=3
    B=2
    B=5
    B=8
    C=77"""
    grammar = ipd.grammar()
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [2,5,8], 'C':77}}, grammar)
    ip=ipd.read_from_string(data)
    self.assertEqual({'ENERGY': {'A' : 3, 'B' : ar([2,5,8]), 'C':77}}, ip.to_dict())
    self.assertEqual('ENERGY A=3 B=2 B=5 B=8 C=77', t(ip.ENERGY.to_string()))
    ipd['ENERGY'].force_order=True
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [2,5,8], 'C':77}}, grammar)
    ip=ipd.read_from_string(data)
    self.assertEqual({'ENERGY': {'A' : 3, 'B' : ar([2,5,8]), 'C':77}}, ip.to_dict())
    self.assertEqual('ENERGY A=3 B=2 B=5 B=8 C=77', t(ip.ENERGY.to_string()))

    ipd["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
    grammar = ipd.grammar()
    data = """ENERGY A=3
    B1=2
    B2=5
    B3=8
    C=77"""
    pp.ParserElement.verbose_stacktrace = True
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [2,5,8], 'C':77}}, grammar)
    self.assertEqual(t(data), t(ip.ENERGY.to_string()) )
    data = """ENERGY A=3
    B1=2
    B3=5
    B2=8
    C=77"""
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [2,8,5], 'C':77}}, grammar)

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('B', gt.Array(int, length=2), is_repeated = True),
        V('C', 3) ]})
    grammar = ipd.grammar()

    data = """ENERGY A=3
    B=2 3
    B=5 9
    B=8 16
    C=77"""
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [ar([2,3]),ar([5,9]),ar([8,16])], 'C':77}}, grammar)
    ip=ipd.read_from_string(data)
    self.assertEqual(t(data), t(ip.ENERGY.to_string()) )

    ipd["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
    grammar = ipd.grammar()
    assertNotValid(data)
    data = """ENERGY A=3
    B1=2 3
    B2=5 9
    B3=8 16
    C=77"""
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : [ar([2,3]),ar([5,9]),ar([8,16])], 'C':77}}, grammar)
    ip=ipd.read_from_string(data)
    self.assertEqual(t(data), t(ip.ENERGY.to_string()) )

#
  def test_sparse_numbered(self):
    input_parameters_def = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        V('NE', gt.SetOf(int, min_length=1)),
        V('Ime', float, 0.0),
      ],
      'SITES' : [
        V('NL', 1)
      ]
    })
    with generate_grammar():
      grammar = input_parameters_def._grammar_of_values()

    assertParse = partial(self.assertParse, grammar=lambda: grammar)
    assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

    assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0")
    assertParse("ENERGY NE=1 Ime=0.5", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : 0.5}})
    input_parameters_def.sections['ENERGY']['Ime'].is_repeated = V.Repeated['DICT']

    with generate_grammar():
      grammar = input_parameters_def._grammar_of_values()
    assertNotValid("ENERGY NE=1 Ime=0.5")
    assertNotValid("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8")
    assertParse("ENERGY NE=1 Ime1=0.4 Ime5=0.8", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : {1:0.4, 5:0.8} }, })

    input_parameters_def.sections['ENERGY']['Ime'].is_repeated = V.Repeated['DEFAULTDICT']
    with generate_grammar():
      grammar = input_parameters_def._grammar_of_values()

    assertParse("ENERGY NE=1 Ime1=0.4 Ime5=0.8", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : {1:0.4, 5:0.8} }, })
    assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0")
    assertParse("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8", { 'ENERGY' : { 'NE' : ar(1), 'Ime' : {'def': 0.5, 1:0.4, 5:0.8} }, })

    ip=input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
    self.assertEqual(0.5, ip.ENERGY.Ime())
    self.assertEqual(0.8, ip.ENERGY.Ime[5])
    self.assertEqual([0.4,0.5,0.5,0.5,0.8], ip.ENERGY.Ime[:])
    self.assertEqual({'def': 0.5, 1: 0.4, 5: 0.8}, ip.ENERGY.Ime.as_dict())
    self.assertEqual({'def': 0.5, 1: 0.4, 5: 0.8}, ip.ENERGY.Ime(all_values=True))
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
    self.assertEqual([0.3,0.2,0.3], ip.ENERGY.Ime[[1,5,9]])
    ip.ENERGY.Ime = { 3: 5., 'def': 3.}
    self.assertEqual([3.,3.,5.], ip.ENERGY.Ime[:])
    ip.ENERGY.Ime[2:5] = [ 2., 3., 8.]
    self.assertEqual([3.,2.,3.,8.], ip.ENERGY.Ime[:])
    ip.ENERGY.Ime[1:3] = 7.
    self.assertEqual([7.,7.,3.,8.], ip.ENERGY.Ime[:])

    ip=input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
    with pytest.raises(DataValidityError):
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

    with pytest.raises(DataValidityError):
      ip.ENERGY.Ime[5] = 'ss'
    ip.ENERGY.Ime.set_dangerous('yy', index=5)
    out=ip.ENERGY.to_string()
    self.assertEqual('ENERGY GRID={3} NE={1} Ime=1.0 Ime1=0.4 Ime5=yy', re.sub(r'\s+',' ', ip.ENERGY.to_string()).strip() )
    with pytest.raises(pp.ParseBaseException):
         input_parameters_def.read_from_file(io.StringIO(out))
    ip=input_parameters_def.read_from_file(io.StringIO(out), allow_dangerous=True)
    self.assertEqual(ip.ENERGY.Ime(all_values=True), {'def': 1.0, 1:0.4, 5:'yy' } )

  def test_gather(self):
    assertParse = self.assertParse
    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('GRID', gt.SetOf(int, length=1), fixed_value=3),
        *gather( V('A', 1),
                V('B', 2) ),
        V('C', 3)
      ]
    })
    assertParse("ENERGY GRID={3} A B=1 2 C=3", {'ENERGY': { 'GRID': ar(3), 'A' : 1, 'B' : 2, 'C' : 3 }}, ipd.grammar())
    out = ipd.read_from_string("ENERGY GRID={3} A B=1 2 C=3")
    # self.assertEqual("ENERGY GRID={3} A B=1 2 C=3 TASK INPUTPARAMETERSDEFINITION ", re.sub(r'[\s\t\n]+',' ', out.to_string()))
    self.assertEqual("ENERGY GRID={3} A B=1 2 C=3 ", re.sub(r'[\s\t\n]+',' ', out.to_string()))

  def test_numpy_array(self):
    assertParse = self.assertParse

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('C', gt.NumpyArray(lines=3), name_in_grammar=False ),
      ]
    })
    assertParse("""ENERGY
    1 1
    2 2
    3 3
    """,{'ENERGY': {'C': ar([1,1,2,2,3,3]).reshape(3,2)}}, ipd.grammar())

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('B', 2),
        V('C', gt.NumpyArray(lines=3), name_in_grammar=False ),
        V('D', gt.NumpyArray(lines=2), name_in_grammar=False),
        V('E', 3),
      ]
    })
    data = """ENERGY A=3 B=2
    1 1
    2 2
    3 3
    4 5 6
    5 9 8
    E=77"""

    g = ipd.grammar()
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : 2, 'C': ar([1.,1,2,2,3,3]).reshape(3,2), 'D': ar([4.,5,6,5,9,8]).reshape((2,3)),'E':77}}, g)
    ipd['ENERGY'].force_order=True
    g = ipd.grammar()
    assertParse(data, {'ENERGY': {'A' : 3, 'B' : 2, 'C': ar([1.,1,2,2,3,3]).reshape(3,2), 'D': ar([4.,5,6,5,9,8]).reshape((2,3)),'E':77}}, g)

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('B', 2),
        V('C', gt.NumpyArray(lines='A'), name_in_grammar=False ),
        V('D', gt.NumpyArray(lines='B'), name_in_grammar=False),
        V('E', 3),
      ]
    })
    assertParse("""ENERGY A=3 B=2
    1 1
    2 2
    3 3
    4 5 6
    5 9 8
    E=77""",{'ENERGY': {'A' : 3, 'B' : 2, 'C': ar([1.,1,2,2,3,3]).reshape(3,2), 'D': ar([4.,5,6,5,9,8]).reshape((2,3)), 'E':77}}, ipd.grammar())

    def test(na, val, str):
        self.assertEqual(na.string(val), str)
        with generate_grammar():
          self.assertEqual(np.asarray(val), na.parse(str))

    test(gt.NumpyArray(item_format='%2.0f'),[[1,2,3],[4,5,6]], " 1  2  3\n 4  5  6")
    test(gt.NumpyArray(indented=2, item_format='%2.0f'), [[1,2,3],[4,5,6]],
                       "   1  2  3\n   4  5  6")
    test(gt.NumpyArray(indented=(8, 2), item_format='%2.0f'),[[1,2,3,4,5,6],[4,5,6,7,8,9]],
                       " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9")
    test(gt.NumpyArray(indented=(8, 2), no_newline_at_end=False, item_format='%2.0f'),[[1,2,3,4,5,6],[4,5,6,7,8,9]],
                       " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9\n")
    test(gt.NumpyArray(line_length=8, item_format='%2.0f', shape=(2,-1)),
                       [[1,2,3,4,5,6],[4,5,6,7,8,9]],
                       " 1  2  3\n 4  5  6\n 4  5  6\n 7  8  9")

  def test_copy(self):
    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('B', 2, is_repeated = True),
        V('C', 3) ]})
    ipd2 = ipd.copy()
    ipd2['ENERGY']['A'].default_value = 5
    assert ipd['ENERGY']['A'].default_value == 1

  def test_length(self):
    t=lambda x:re.sub(r'\s+',' ',x).strip()

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('NK', Length('K1', 'K2')),
        V('K1', gt.Array(int), is_optional=True),
        V('K2', gt.Array(int), is_optional=True)
      ]
    })
    ip = ipd.create_object()
    with pytest.warns(DataValidityError):
        ip.ENERGY.K1 = [1,2,3]
    ip.ENERGY.K2 = [2,2,3]
    assert t(ip.to_string()) == 'ENERGY A=1 NK=3 K1=1 2 3 K2=2 2 3'
    ip2 = ipd.read_from_string(ip.to_string())
    assert ip2.ENERGY.NK() == 3

    with pytest.raises(DataValidityError):
        ip2.ENERGY.NK=7

    with pytest.warns(DataValidityError):
        ip.ENERGY.K2 = [2,2]
    with pytest.raises(DataValidityError):
        ip.to_string()

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('NK', Length('K1', 'K2', default_values=(1,2))),
        V('K1', gt.Array(int), is_optional=True),
        V('K2', gt.Array(int), is_optional=True)
      ]
    })
    ip = ipd.create_object()
    ip.ENERGY.NK = 3
    assert (ip.ENERGY.K1() == [1,1,1]).all()
    assert (ip.ENERGY.K2() == [2,2,2]).all()
    ip.ENERGY.NK = 5
    assert (ip.ENERGY.K1() == [1,1,1,1,1]).all()
    ip.ENERGY.NK = 2
    assert (ip.ENERGY.K1() == [1,1]).all()
    ip.ENERGY.NK = 2

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        V('NK', Length('K1', default_values=[1,2,3])),
        V('K1', gt.Array(int, length=3), is_repeated='REPEATED', is_optional=True),
      ]
    })
    ip = ipd.create_object()
    ip.ENERGY.NK = 2
    assert (ip.ENERGY.K1() == [[1,2,3],[1,2,3]]).all()
    ip.ENERGY.K1[1]= [3,3,3]
    ip.ENERGY.NK = 3
    assert (ip.ENERGY.K1() == [[1,2,3],[3,3,3],[1,2,3]]).all()

#
  def test_switch(self):
    assertParse = self.assertParse
    assertNotValid = self.assertNotValid
    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        *switch('A', {
          1: [ V('B', 2) ],
          2: V('C', 3) })
       ]})

    ipd.custom_class=None
    ipd['ENERGY'].custom_class=None
    ipd['ENERGY'].force_order=True
    grammar = ipd.grammar()
    assertParse('ENERGY A=2 C=4', {'ENERGY': {'A' : 2, 'C' : 4}}, grammar)
    assertParse('ENERGY A=2 C=4', {'ENERGY': {'A' : 2, 'C' : 4}}, grammar)
    out=ipd.read_from_string('ENERGY A=2 C=4')
    self.assertEqual(out['ENERGY'].to_string(), 'ENERGY\n\tA=2\n\tC=4\n')
    self.assertEqual(out['ENERGY'].to_string(), 'ENERGY\n\tA=2\n\tC=4\n')
    assertParse('ENERGY A=1 B=2', {'ENERGY': {'A' : 1, 'B' : 2}}, grammar)
    out=ipd.read_from_string('ENERGY A=1 B=4')
    self.assertEqual(out['ENERGY'].to_string(), 'ENERGY\n\tA=1\n\tB=4\n')

    assertNotValid('ENERGY A=2 B=2', grammar)
    assertNotValid('ENERGY A=1 C=1', grammar)
    assertNotValid('ENERGY A=1 B=2 C=1', grammar)
    assertNotValid('ENERGY A=2 B=2 C=1', grammar)
    assertNotValid('ENERGY A=3 C=3', grammar)

    ipd = cd.InputParametersDefinition.definition_from_dict({
      'ENERGY' : [
        V('A', 1),
        *switch('A', {
          1: [ V('B', 2), V('D', 3), V('C',1) ],
          2: [ V('C', 3), V('E', 1), V('B',2)] })
       ]})
    ipd['ENERGY'].force_order = True
    grammar = ipd.grammar()
    assertParse('ENERGY A=2 C=4 E=2 B=6', {'ENERGY': {'A':2, 'C':4, 'E':2, 'B':6 }}, grammar)
    assertParse('ENERGY A=1 B=5 D=3 C=5', {'ENERGY': {'A':1, 'B':5, 'D':3, 'C':5 }}, grammar)

    grammar = ipd.grammar()
    assertNotValid('ENERGY A=2 B=2', grammar)
    assertNotValid('ENERGY A=1 B=2', grammar)
    assertNotValid('ENERGY A=1 B=2 D=3', grammar)
    assertNotValid('ENERGY A=1 B=5 D=3 C=5 B=6', grammar)
    assertNotValid('ENERGY A=1 B=5 D=3 C=5 E=6', grammar)
    assertNotValid('ENERGY A=2 B=5 C=4 E=3', grammar)
    out=ipd.read_from_string('ENERGY A=2 C=4 E=2 B=6')
    self.assertEqual(out['ENERGY'].to_dict(), {'A': 2, 'C': 4, 'E':2, 'B':6})
    self.assertEqual(out['ENERGY'].to_string(), 'ENERGY\n\tA=2\n\tC=4\n\tE=2\n\tB=6\n')
