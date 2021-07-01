if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..potential_definitions import PotentialDefinition
from ..custom_potential_section import SectionString
from ...common.grammar import generate_grammar
import re
from ...common.grammar import delimitedList

class TestCustomSection(TestCase):

  def assertNone(self, val):
      return self.assertEqual(None, val)

  def test_custom_section(self):

    with generate_grammar():
      self.assertTrue(re.compile(SectionString.delimiter_pattern).fullmatch('\n****************\n'))
      self.assertTrue(re.compile(SectionString.delimiter_pattern).fullmatch('    \n****************\n'))
      self.assertTrue(re.compile(SectionString.delimiter_pattern).fullmatch(' \n \n   \n****************\n'))
      self.assertTrue(re.compile(SectionString.delimiter_pattern).fullmatch('\n****************\n   \n    \n'))
      self.assertNone(re.compile(SectionString.delimiter_pattern).fullmatch('a\n****************\n'))
      self.assertNone(re.compile(SectionString.delimiter_pattern).fullmatch('\n****************\na'))
      self.assertNone(re.compile(SectionString.delimiter_pattern).fullmatch('\n******* *********\na'))

      cmg = PotentialDefinition.custom_member_grammar()

    sec = """HOST MADELUNG POTENTIAL
        IQ      VLMMAD
NLMTOP-POT         4
         1   1 -2.34816691089662E-01
         1   2  0.00000000000000E+00to je konec"""

    cmg.parseString(sec, True)
    self.assertTrue(cmg.parseString(sec + '\n')[0][1].endswith('to je konec'))

    with generate_grammar():
      cmgs = delimitedList(cmg, SectionString.grammar_of_delimiter())
    out = cmgs.parseString(sec + "\n************************\n" + sec, True)
    self.assertEqual(2, len(out))


    out = cmgs.parseString(sec +
                          "\n************************\n" +
                          sec + "     \n****************************       \n   \n" +
                          sec + sec,
                            True)
    self.assertEqual(3, len(out))
