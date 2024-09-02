"""
Wrapper for spglib for computing symmetry of a primitive cell
"""
from __future__ import annotations

from typing import Optional, List, Dict, Tuple
import spglib
from ase import Atoms
from ase.spacegroup import Spacegroup
from ..common.unique_values import UniqueValuesMapping
from ..sprkkr import sprkkr_atoms
import numpy as np


class SpacegroupInfo:
    """ Class, that carry information about spacegroup and symmetry of a structure """

    def __init__(self, atoms: sprkkr_atoms.SPRKKRAtoms,
                       spacegroup: Optional[Spacegroup],
                       dataset: Optional[Dict]=None,
                       equivalent_sites: Optional[UniqueValuesMapping]=None):
        """
        Parameters
        ----------
        atoms
          ASE atoms object desribed by the SpacegroupInfo

        spacegroup
          ASE spacegroup object. If None, there is no symmetry in the crystal.

        dataset
          SPGLib dataset - a dictionary, containing many informations about the symmetry and
          spacegroup of the structure. If it is not provided, and ``spacegroup`` is not ``None``
          the spacegroup is recomputed to obtain the dataset, if the dataset is requested.

        equivalent_sites
          Object, that provides equivalence classes for the atoms. If it is not provided,
          and ``spacegroup`` is not ``None``, the spacegroup is recomputed to obtain the equivalence
          equivalent_sites, if the equivalent_sites is requested.
        """

        atoms = sprkkr_atoms.SPRKKRAtoms.promote_ase_atoms(atoms)
        self.atoms = atoms
        self.spacegroup = spacegroup
        self._dataset = dataset
        self._equivalent_sites = equivalent_sites

    def __repr__(self):
        return str(self.spacegroup) or 'No spacegroup'

    def __str__(self):
        return self.__repr__()

    def number(self) -> Optional[int]:
        """
        Returns
        -------
        spacegroup
          Spacegroup number or None, if there is no spacegroup.
        """
        return self.spacegroup.no if self.spacegroup else None

    @property
    def dataset(self)->Optional[Dict]:
        if self._dataset is None and self.spacegroup:
           self.recompute(allow_change=False)

        return self._dataset

    @property
    def equivalent_sites(self)->UniqueValuesMapping:
        """
          Object, that provide equivalence classes for the atoms. If it is not provided,
          and ``spacegroup`` is not ``None``, the spacegroup is recomputed to obtain the equivalent_sites,
          if the equivalent_sites is requested.
        """
        if self._equivalent_sites is None:
           if self._dataset:
              return UniqueValuesMapping(self._dataset['equivalent_atoms'])
           elif self.spacegroup:
              self.recompute(allow_change=False)
           else:
              self._equivalent_sites = UniqueValuesMapping(np.arange(len(self.atoms)))
        return self._equivalent_sites

    def recompute(self, allow_change:bool=True):
        """ Recompute spacegroup, dataset and equivalent_sites.

        allow_change
          If False, the resulting spacegroup have to be the same as the old one.
          Used, if the spacegroup info have been constructed using a spacegroup
          without providing spglib dataset.
        """

        if not self.atoms.symmetry:
           self.spacegroup = None
           self._dataset = None
           self._equivalent_sites = UniqueValuesMapping(np.arange(len(self.atoms)))

        sg, dataset, equivalent_sites = self.compute_spacegroup_info(self.atoms)
        if not allow_change and (
             bool(sg) != bool(self.spacegroup) or
             (sg and self.spacegroup.no != sg.no)):
                 raise "The stored spacegroup do not correspond to the actual one"

        self.spacegroup = sg
        self._dataset = dataset
        self._equivalent_sites = equivalent_sites

    @staticmethod
    def compute_spacegroup_info(atoms:Atoms,
                       atomic_numbers : List=None,
                       consider_old   : bool=False,
                       precision      : float=1e-5,
                       angular_precision:float=1.0,
                       ) -> Tuple[Spacegroup, Dict, UniqueValuesMapping]:

       """ Return the values needed to create (or update) :class:`Spacegroup`
           that suits to the atoms' cell structure.

           For the parameters description, see the ':meth:`from_atoms` method.

           Returns
           -------
           spacegroup
              ASE Spacegroup object, describing the spacegroup

           dataset
              Dataset, describing symmetry- and spacegroup-related information
              about the structure, see the `Spglib documentation <https://spglib.github.io/spglib/dataset.html#spglib-dataset>`
              for the documentation of the members of the `dataset and the `documentation of the C-to-python Spglib binding
              <https://spglib.github.io/spglib/python-spglib.html>` for the names of the members (which slightly differ from
              the names in the C structure)

           equivalent_sites
              The map of equivalent atoms.
       """

       sg_dataset = None
       sprkkr_atoms.SPRKKRAtoms.promote_ase_atoms(atoms)
       if atoms.symmetry:
           equivalent_sites = possibly_equivalent_sites(atoms, atomic_numbers, consider_old)
           spositions = atoms.get_scaled_positions()
           sg_dataset = spglib.get_symmetry_dataset((atoms.get_cell(),
                             spositions,
                             equivalent_sites.mapping),
                             symprec = precision,
                             angle_tolerance = angular_precision)
       if sg_dataset is None:
          return None, None, UniqueValuesMapping(np.arange(len(atoms)))

       spacegroup = Spacegroup(sg_dataset['number'])
       equivalent_sites = equivalent_sites.merge(sg_dataset['equivalent_atoms'])
       equivalent_sites.normalize(start_from=0)
       return spacegroup, sg_dataset, equivalent_sites

    @staticmethod
    def from_atoms(atoms:Atoms,
                       atomic_numbers : List=None,
                       consider_old   : bool=False,
                       precision      : float=1e-5,
                       angular_precision:float=1.0,
                       ) -> SpacegroupInfo:
        """ Create :class:`SpacegroupInfo` for given ASE :class:`Atoms` object.

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
        """

        out = SpacegroupInfo.compute_spacegroup_info(**locals())
        return SpacegroupInfo(atoms, *out)


def possibly_equivalent_sites(atoms: Atoms,
                       atomic_numbers : List=None,
                       consider_old   : bool=False) -> UniqueValuesMapping:
    """ Return the object that describe equivalence classes of the atoms
    from the given Atoms object (i.e. in one equivalence class will be the
    atoms with the same atomic number, occupation, etc.)
    """
    if atomic_numbers is not None:
        equivalent_sites = UniqueValuesMapping(atomic_numbers)
    else:
        equivalent_sites = UniqueValuesMapping(atoms.get_atomic_numbers())

    if consider_old:
        try:
           equivalent_sites = equivalent_sites.merge([i.site_type if i else None for i in atoms.get_array(atoms.sites_array_name)])
        except KeyError:
           pass
    else:
        occupation = atoms.info.get('occupancy', {})
        if occupation:
            def gen_occ():
              for i in range(len(equivalent_sites)):
                  val = occupation.get(i, None) or \
                        occupation.get(str(i), None)
                  if val is None:
                      yield val
                  else:
                      yield tuple((k, val[k]) for k in val)
            equivalent_sites = equivalent_sites.merge(gen_occ())
    equivalent_sites.normalize(start_from=0)
    return equivalent_sites
