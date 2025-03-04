from ase.build import bulk
import pytest

if __package__:
   from .init_tests import TestCase, patch_package
else:
   from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..spglib import spglib_dataset  # NOQA: E402
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms  # NOQA: E402
from ...common.unique_values import UniqueValuesMapping  # NOQA: E402


class TestSpgilib(TestCase):

  def test_spglib_dataset(self):

       def ssert(mapping, out):
           self.assertTrue(UniqueValuesMapping(out.equivalent_atoms).is_equivalent_to(mapping))

       atoms = bulk('NaCl', "rocksalt", a=5.64)
       dataset = spglib_dataset(atoms)
       ssert([0,1], dataset)
       ssert([0,0], spglib_dataset(atoms, atomic_numbers=[4,4]))

       with pytest.raises(AttributeError):
            dataset.non_existent_value

       atoms = bulk('NaNa', "rocksalt", a=5.64)
       ssert([0,0], spglib_dataset(atoms))
       ssert([0,1], spglib_dataset(atoms, atomic_numbers=[4,7]))

       atoms = bulk('NaCl', "rocksalt", a=5.64)
       SPRKKRAtoms.promote_ase_atoms(atoms)
       atoms.sites
       ssert([0,1], spglib_dataset(atoms))
       atoms.sites[1].occupation = {'Na':1.0}
       ssert([0,1], spglib_dataset(atoms, consider_old=True))
       # repeat the test to be sure
       ssert([0,1], spglib_dataset(atoms, consider_old=True))
       ssert([0,0], spglib_dataset(atoms, consider_old=False))
