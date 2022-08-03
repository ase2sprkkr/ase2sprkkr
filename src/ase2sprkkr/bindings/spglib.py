"""
Wrapper for spglib for computing symmetry of a primitive cell
"""
from ase import Atoms
from typing import List, Dict

import spglib
from ase.spacegroup import Spacegroup
from ..common.unique_values import UniqueValuesMapping

def possibly_equivalent_sites(atoms: Atoms,
                       atomic_numbers : List=None,
                       consider_old   : bool=False) -> UniqueValuesMapping:
    """ Return the object that describe equivalence classes of the atoms
    from the given Atoms object (i.e. in one equivalence class will be the
    atoms with the same atomic number, occupation, etc.)
    """
    if atomic_numbers is not None:
        mapping = UniqueValuesMapping(atomic_numbers)
    else:
        mapping = UniqueValuesMapping(atoms.get_atomic_numbers())
    if consider_old:
        try:
           mapping = mapping.merge(atoms.get_array(atoms.sites_array_name))
        except KeyError:
           pass
    else:
        occupation = atoms.info.get('occupancy', {})
        if occupation:
            def gen_occ():
              for i in range(len(mapping)):
                  val = occupation.get(i, None)
                  if val is None:
                      yield val
                  else:
                      yield tuple((k, val[k]) for k in val)
            mapping = mapping.merge(gen_occ())
    return mapping


def compute_spacegroup(atoms:Atoms,
                       atomic_numbers : List=None,
                       consider_old   : bool=False,
                       precision      : float=1e-5,
                       angular_precision:float=1.0) -> (Spacegroup, Dict):

       """ Return the spacegroup that suits to the atoms' cell structure.

           Parameters
           ----------
           atoms
            ASE atoms structure

           atomic_numbers
            The atomic numbers can be overriden, e.g. to force (partialy) break the symmetry.
            The atomic numbers have not to be the real ones, they can be (and are treated) just 
            the labels of equivalence classes.

           consider_old

           precision
            The tolerated unprecision

           angular_precision
            The tolerated unprecision in angular coordinate

            Returns
            -------
            spacegroup
              ASE spacegroup object

            dataset
              Spglib symmetry dataset
       """

       mapping = possibly_equivalent_sites(atoms, atomic_numbers, consider_old)
       spositions = atoms.get_scaled_positions()
       sg_dataset = spglib.get_symmetry_dataset((atoms.get_cell(),
                             spositions,
                             mapping.mapping),
                             symprec = precision,
                             angle_tolerance = angular_precision)
       if sg_dataset is None:
           return None
       spacegroup = Spacegroup(sg_dataset['number'])
       spacegroup.dataset = sg_dataset
       tags = spacegroup.tag_sites(spositions)
       mapping = mapping.merge(tags)
       spacegroup.sites_equivalence_mapping = mapping
       return spacegroup

def equivalent_sites_for_spacegroup(atoms, spacegroup, atomic_numbers=None, consider_old=False):
    """
    Return class that describe equivalent classes for a given spacegroup. If the spacegroup is not
    created by :func:`compute_spacegroup`, the precomputed equivalent classes are taken from it.
    Otherwise, the additional parameters are used to determine the  equivalent atoms.
    """
    try:
      mapping = spacegroup.sites_equivalence_mapping
    except AttributeError:
      mapping = possible_equivalent_atoms(atoms, atomic_numbers=None, consider_old=consider_old)
    spositions = atoms.get_scaled_positions()
    return mapping.merge(spacegroup.tag_sites(spositions))

def spacegroup_dataset(atoms):
    if not atoms.symmetry:
       return None
    sg = atoms.info.get('spacegroup', None)
    if sg:
       if not hasattr(sg, 'dataset'):
          #spacegroup has not been created by ASE2SPRKKR, recompute it to obtain the dataset
          sg_new = compute_spacegroup(atoms)
          if sg.no != sg_new.no:
             raise "The stored spacegroup do not correspond to the actual one"
          atoms.info['spacegroup'] = sg = sg_new
    else:
       atoms.compute_sites_symmetry()
       if not atoms.symmetry:
          return None

    return sg.dataset
