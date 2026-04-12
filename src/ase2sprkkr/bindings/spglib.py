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
from ..sprkkr.occupations import Occupation

if hasattr(spglib, "set_error_handling"):
    spglib.set_error_handling(False)

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
      If True, tag according to the spacegroup_kinds (which foor SPRKKRAtoms
      reflects the site_types.
      If False, tag just according to the atomic numbers and, for SPRKKRAtoms,
      occupancy (thus, possibly merge the old differnt site_types into one).
    """
    def mapp(x):
        if return_mapping or consider_old:
            return UniqueValuesMapping(x)
        else:
            return x

    def unmap(x):
        return x if return_mapping else x.mapping

    if hasattr(atoms, 'are_sites_inited') and atoms.are_sites_inited():
        if consider_old and np.sum(atoms.arrays['sprkkr_sites_data'] == 0) == 0:
            out = [ i.site_type for i in atoms.sites ]
        else:
            def occ(i, site):
                return site.occupation if site else Occupation.for_atom_from_ase_atoms(atoms, i)
            out = ( occ(i, site).as_dict for i, site in enumerate(atoms.sites) )
            out = [str(sorted(i.items())) for i in out]
    else:
        out = atoms.get_atomic_numbers()

    if consider_old and 'spacegroup_kinds' in atoms.arrays:
        return unmap(UniqueValuesMapping(out).
                  merge(atoms.arrays['spacegroup_kinds']))
    else:
        return mapp(out)




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

    def dataset(equivalent_sites):
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

    sg_dataset = dataset(equivalent_sites)
    if sg_dataset and consider_old and \
        hasattr(atoms, 'are_sites_inited') and atoms.are_sites_inited() and \
        np.sum(atoms.arrays['sprkkr_sites_data'] == 0) != 0:
            equivalent_sites = sg_dataset.equivalent_atoms
            recompute = False
            types = {}
            counter = max(equivalent_sites)
            for i,(site,typ) in enumerate(zip(atoms.arrays['sprkkr_sites_data'], equivalent_sites)):
                if site == 0:
                    continue
                stype = site.site_type
                if not typ in types:
                    types[typ] = { stype : typ }
                    out = typ
                else:
                    out = types[typ].get(stype, None)
                if out is not None:
                    equivalent_sites[i] = out
                else:
                    counter+=1
                    types[typ][stype] = counter
                    equivalent_sites[i] = counter
            sg_dataset = dataset(equivalent_sites)

    return sg_dataset
