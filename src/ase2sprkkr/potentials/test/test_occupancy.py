import ase
import ase2sprkkr
from ase.build import bulk


if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
   from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
   from ...sprkkr.calculator import SPRKKR


class TestOccupancy(TestCase):

    def test(self):

        print('ASE version:', ase.__version__)
        print('ASE2SPRKKR version:', ase2sprkkr.__version__)

        Cu = bulk('Cu',a=3.6)

        atoms = SPRKKRAtoms.promote_ase_atoms(Cu)
        atoms.sites[0].occupation.set({'Cu':0.75, 'Ni':0.25 })

        scf_opts = {
           'SITES.NL':3,
           'TAU.BZINT':'POINTS',
           'TAU.NKTAB3D':500,
           'SCF.VXC':'VWN',
           'SCF.NITER':200,
           'SCF.MIX':0.1,
           'SCF.MIXOP':0.1,
           'SCF.TOL':1E-5,
        }

        calculator=SPRKKR(atoms=atoms)
        calculator.save_input(input_file='Cu.inp', potential_file='Cu.pot', options=scf_opts)
