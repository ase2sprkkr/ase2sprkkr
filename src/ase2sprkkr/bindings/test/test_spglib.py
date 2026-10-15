from ase.build import bulk
import pytest

from ..spglib import spglib_dataset
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ...common.unique_values import UniqueValuesMapping


def assert_equivalent(mapping, dataset):
    assert UniqueValuesMapping(dataset.equivalent_atoms).is_equivalent_to(mapping)


def test_spglib_dataset():
    atoms = bulk("NaCl", "rocksalt", a=5.64)
    dataset = spglib_dataset(atoms)
    assert_equivalent([0, 1], dataset)
    assert_equivalent([0, 0], spglib_dataset(atoms, atomic_numbers=[4, 4]))

    with pytest.raises(AttributeError):
        dataset.non_existent_value

    atoms = bulk("NaNa", "rocksalt", a=5.64)
    assert_equivalent([0, 0], spglib_dataset(atoms))
    assert_equivalent([0, 1], spglib_dataset(atoms, atomic_numbers=[4, 7]))

    atoms = bulk("NaCl", "rocksalt", a=5.64)
    SPRKKRAtoms.promote_ase_atoms(atoms)
    atoms.sites
    assert_equivalent([0, 1], spglib_dataset(atoms))
    atoms.sites[1].occupation = {"Na": 1.0}
    assert_equivalent([0, 1], spglib_dataset(atoms, consider_old=True))
    assert_equivalent([0, 1], spglib_dataset(atoms, consider_old=True))
    assert_equivalent([0, 0], spglib_dataset(atoms, consider_old=False))
