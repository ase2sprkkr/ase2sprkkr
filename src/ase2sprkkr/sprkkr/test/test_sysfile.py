if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ase.build import bulk
from ase.spacegroup import crystal
from ..sysfile import sysfile_content, write_sysfile
import os, tempfile

class SysfileTest(TestCase):

    def test(self):
        a = 4.6
        c = 2.95
        atoms = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                    spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
        #atoms = bulk('Ag')
        x=sysfile_content(atoms)
        y="""
system data-file created by python ase2sprkkr
<unknown>
xband-version
5.0
dimension
3D
Bravais lattice
8 tetragonal primitive 4/mmm D_4hspace group number (ITXC and AP)
  136  384
structure type
UNKNOWN
lattice parameter A  [a.u.]
    8.692740178850
ratio of lattice parameters  b/a  c/a
    1.000000000000    0.641304347826
lattice parameters  a b c  [a.u.]
    8.692740178850    8.692740178850    5.574692071219
lattice angles  alpha beta gamma  [deg]
   90.000000000000   90.000000000000   90.000000000000
primitive vectors     (cart. coord.) [A]
    1.000000000000    0.000000000000    0.000000000000
    0.000000000000    1.000000000000    0.000000000000
    0.000000000000    0.000000000000    1.000000000000
number of sites NQ
  6
 IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ
  1   1    0.000000000000    0.000000000000    0.000000000000      0.000000000000   4    1   1
  2   1    2.300000000000    2.300000000000    1.475000000000      0.000000000000   4    1   1
  3   2    1.380000000000    1.380000000000    0.000000000000      0.000000000000   4    1   2
  4   2    3.220000000000    3.220000000000    0.000000000000      0.000000000000   4    1   2
  5   2    0.920000000000    3.680000000000    1.475000000000      0.000000000000   4    1   2
  6   2    3.680000000000    0.920000000000    1.475000000000      0.000000000000   4    1   2
number of sites classes NCL
  2
ICL WYCK NQCL IQECL (equivalent sites)
  1   -    2  1  2
  2   -    4  3  4  5  6
number of atom types NT
  2
 IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)
  1  22        Ti    2 1.000  1  2
  2   8         O    4 1.000  3  4  5  6
        """
        sanitize = lambda x: "\n".join(map(lambda x: x.strip(), x.strip().split('\n')))
        self.assertEqual(sanitize(x),sanitize(y))
        if os.name != 'nt':
           write_sysfile(atoms, tempfile.NamedTemporaryFile().name)
