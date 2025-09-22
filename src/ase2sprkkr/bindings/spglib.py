"""
Wrapper for spglib for computing symmetry of a primitive cell
"""
from __future__ import annotations
from typing import List, Union, Optional
from ..sprkkr.atoms_region import AtomsRegion
import spglib
import numpy as np
from ase import Atoms
from ..common.unique_values import UniqueValuesMapping


def spglib_dataset_wrapper(dataset):
    """ Backward compatibility """
    if dataset is None or \
       hasattr(dataset, 'equivalent_atoms'):
         return dataset

    class Convertor:

         def __getattr__(self, name):
             try:
                 return dataset[name]
             except KeyError as ke:
                 raise AttributeError(ke)

    return Convertor()


def tag_sites(atoms:Atoms, consider_old:bool, return_mapping:False):
    """
    Return array of class-equivalence-numbers of sites of a
    given Atoms object

    atoms
      The atoms to return the equivalence classes

    consider_old
      If True, tag according to the spacegroup_kinds (which for SPRKKRAtoms
      reflects the site_types.
      If False, tag just according to the atomic numbers and, for SPRKKRAtoms,
      occupancy (thus, possibly merge the old differnt site_types into one).
    """
    def mapp(x):
        if return_mapping:
            return UniqueValuesMapping(x)
        else:
           return x

    def unmap(x):
        return x if return_mapping else x.mapping

    if hasattr(atoms, 'are_sites_inited') and atoms.are_sites_inited() and \
                  np.sum(atoms.arrays['sprkkr_sites_data'] == 0) == 0:
                  # if SPRKKR Atoms is added to Atoms, it seems that has the sites inited
        if consider_old:
            return mapp(atoms.arrays['spacegroup_kinds'])
        out = ( i.site_type.occupation.as_dict for i in atoms.sites )
        return mapp(np.unique([str(sorted(i.items())) for i in out], return_inverse=True)[1])

    equivalent_sites = atoms.get_atomic_numbers()
    if not ('spacegroup_kinds' in atoms.arrays and consider_old):
        return mapp(equivalent_sites)
    return unmap(UniqueValuesMapping(equivalent_sites).
        merge(atoms.arrays['spacegroup_kinds']))


def spglib_dataset(atoms: "Union[Atoms,AtomsRegion]",
                       atomic_numbers : Optional[List]=None,
                       consider_old: bool = True,
                       precision      : float=1e-5,
                       angular_precision:float=1.0,
                       add: Optional[List] = None
                       ) -> np.ndarray:
    """ Return the object that describe equivalence classes of the atoms
    from the given Atoms object (i.e. in one equivalence class will be the
    atoms with the same atomic number, occupation, etc.)
    """
    if atomic_numbers is not None:
        dtype = getattr(atomic_numbers, "dtype", None)
        if np.issubdtype(dtype, np.integer):
            equivalent_sites = atomic_numbers
        else:
            equivalent_sites = UniqueValuesMapping(atomic_numbers)
    else:
        equivalent_sites = tag_sites(atoms, consider_old, return_mapping=True)

    if add is not None:
        if not hasattr(equivalent_sites, 'mapping'):
            equivalent_sites = UniqueValuesMapping(equivalent_sites)
        equivalent_sites = equivalent_sites.merge(add).normalize(start_from=0).mapping
    else:
        if hasattr(equivalent_sites, 'mapping'):
            equivalent_sites=equivalent_sites.normalize(start_from=0).mapping

    sg_dataset = spglib.get_symmetry_dataset((atoms.get_cell(),
                       atoms.get_scaled_positions(),
                       equivalent_sites),
                       symprec = precision,
                       angle_tolerance = angular_precision)
    if sg_dataset:
       dataset = spglib_dataset_wrapper(sg_dataset)
    else:
       dataset = False
    return dataset
