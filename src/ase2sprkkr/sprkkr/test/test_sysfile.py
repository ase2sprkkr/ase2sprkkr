from ase.cell import Cell
from ase.spacegroup import crystal

from ...physics.lattice_data import ap_number_from_spacegroup
from ..sysfile import sysfile_content, write_sysfile


def test_ap_number_from_spacegroup():
    assert ap_number_from_spacegroup(3, 3) == 4
    assert ap_number_from_spacegroup(3, 4) == 3
    assert ap_number_from_spacegroup(3, 5) == 3

    assert ap_number_from_spacegroup(5, 9) == 10
    assert ap_number_from_spacegroup(5, 12) == 8
    assert ap_number_from_spacegroup(5, 16) == 7

    assert ap_number_from_spacegroup(7, 21) == 18
    assert ap_number_from_spacegroup(7, 24) == 16
    assert ap_number_from_spacegroup(7, 29) == 15

    assert ap_number_from_spacegroup(9, 39) == 33
    assert ap_number_from_spacegroup(9, 45) == 29
    assert ap_number_from_spacegroup(9, 54) == 27

    assert ap_number_from_spacegroup(136, 419) == 384

    assert ap_number_from_spacegroup(143, cell=Cell.fromcellpar([3, 3, 5, 90, 90, 120])) == 395
    assert ap_number_from_spacegroup(143, cell=Cell.fromcellpar([3, 3, 3, 75, 75, 75])) == 396

    assert ap_number_from_spacegroup(146, 433) == 402
    assert ap_number_from_spacegroup(146, 434) == 401

    assert ap_number_from_spacegroup(160, 450) == 420
    assert ap_number_from_spacegroup(160, 451) == 419

    assert ap_number_from_spacegroup(166, 458) == 428
    assert ap_number_from_spacegroup(166, 459) == 427

    assert ap_number_from_spacegroup(136) == 384


def test_sysfile_content(tmp_path):
    a = 4.6
    c = 2.95
    atoms = crystal(
        ["Ti", "O"],
        basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
        spacegroup=136,
        cellpar=[a, a, c, 90, 90, 90],
    )
    actual = sysfile_content(atoms)
    expected = """
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

    def sanitize(value):
        return "\n".join(line.strip() for line in value.strip().split("\n"))

    assert sanitize(actual) == sanitize(expected)

    output = tmp_path / "system.sys"
    write_sysfile(atoms, str(output))
    assert output.exists()
