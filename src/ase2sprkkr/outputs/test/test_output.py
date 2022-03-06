if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

import os
from ..readers.scf import ScfOutputReader

class TestOutput(TestCase):

  def test_output(self):
      ScfOutputReader.atoms_conf_type.parse(
"""  33 E= 0.6083 0.0000          IT=   1  Li_1
         DOS      NOS     P_spin   m_spin    P_orb    m_orb    B_val      B_core
  s    0.4387   0.0296    0.0000   0.0000   0.00000  0.00000    0.00 s      0.00
  p    1.2579   0.0962    0.0000   0.0000   0.00000  0.00000    0.00 ns     0.00
  d    0.6886   0.0926    0.0000   0.0000   0.00000  0.00000    0.00 cor    0.00
  f    0.3427   0.0476    0.0000   0.0000   0.00000  0.00000    0.00
 sum   2.7279   0.2660    0.0000   0.0000   0.00000  0.00000    0.00 v+c    0.00
 E_band         0.11559127 [Ry]
dipole moment   1      0.0000000000000000      0.0000000000000000      0.0000000000000000""")
