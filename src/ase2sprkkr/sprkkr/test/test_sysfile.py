from ase.spacegroup import crystal
from ase.cell import Cell
import os
import tempfile

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..sysfile import sysfile_content, write_sysfile  # NOQA: E402
from ...physics.lattice_data import ap_number_from_spacegroup  # NOQA: E402


class TestSysfile(TestCase):
    def test_ap_number_from_spacegroup(self):
        self.assertEqual(ap_number_from_spacegroup(3, 3), 4)
        self.assertEqual(ap_number_from_spacegroup(3, 4), 3)
        self.assertEqual(ap_number_from_spacegroup(3, 5), 3)

        self.assertEqual(ap_number_from_spacegroup(5, 9), 10)
        self.assertEqual(ap_number_from_spacegroup(5, 12), 8)
        self.assertEqual(ap_number_from_spacegroup(5, 16), 7)

        self.assertEqual(ap_number_from_spacegroup(7, 21), 18)
        self.assertEqual(ap_number_from_spacegroup(7, 24), 16)
        self.assertEqual(ap_number_from_spacegroup(7, 29), 15)

        self.assertEqual(ap_number_from_spacegroup(9, 39), 33)
        self.assertEqual(ap_number_from_spacegroup(9, 45), 29)
        self.assertEqual(ap_number_from_spacegroup(9, 54), 27)

        self.assertEqual(ap_number_from_spacegroup(136, 419), 384)

        self.assertEqual(ap_number_from_spacegroup(143, cell=Cell.fromcellpar([3, 3, 5, 90, 90, 120])), 395)
        self.assertEqual(ap_number_from_spacegroup(143, cell=Cell.fromcellpar([3, 3, 3, 75, 75, 75])), 396)

        self.assertEqual(ap_number_from_spacegroup(146, 433), 402)
        self.assertEqual(ap_number_from_spacegroup(146, 434), 401)

        self.assertEqual(ap_number_from_spacegroup(160, 450), 420)
        self.assertEqual(ap_number_from_spacegroup(160, 451), 419)

        self.assertEqual(ap_number_from_spacegroup(166, 458), 428)
        self.assertEqual(ap_number_from_spacegroup(166, 459), 427)

        self.assertEqual(ap_number_from_spacegroup(136), 384)

    def test(self):
        a = 4.6
        c = 2.95
        atoms = crystal(["Ti", "O"], basis=[(0, 0, 0), (0.3, 0.3, 0.0)], spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
        # atoms = bulk('Ag')
        x = sysfile_content(atoms)
        y = """
system data-file created by python ase2sprkkr
<unknown>
xband-version
5.0
dimension
3D
Bravais lattice
8 tetragonal primitive 4/mmm D_4h
space group number (ITXC and AP)
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

        def sanitize(x):
            return "\n".join(map(lambda x: x.strip(), x.strip().split("\n")))

        self.assertEqual(sanitize(x), sanitize(y))
        if os.name != "nt":
            write_sysfile(atoms, tempfile.NamedTemporaryFile().name)
