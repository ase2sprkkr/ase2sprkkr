"""
Wrapper for spglib for computing symmetry of a primitive cell
"""
from __future__ import annotations
from typing import Optional, Dict
import numpy as np
import copy
from contextlib import contextmanager

from . import sprkkr_atoms
from ..bindings.spglib import spglib_dataset
from .sites import Site, SiteType


def site_type_copier(atoms):
    """ Return a function, that always copy the site_type """
    def copy(site):
        return site.site_type.copy(atoms=atoms)

    return copy


def used_site_type_copier(atoms):
    """ Return a function, that copy the site_type, if it is requested more than once """
    used = set()

    def copy(site):
        site_type = site.site_type
        if site_type in used or site_type.atoms is not atoms:
            site_type = site_type.copy(atoms=atoms)
        used.add(site_type)
        return site_type

    return copy


class SpacegroupInfo:
    """ Class, that carry information about spacegroup and symmetry of a structure """

    def __init__(self, atoms: sprkkr_atoms.SPRKKRAtoms,
                 symmetry = True,
                 dataset: Optional[Dict]=None):
        """
        Parameters
        ----------
        atoms
          ASE atoms object desribed by the SpacegroupInfo

        symmetry
          Whether the by-symmetry-equivalent sites will share properties

        dataset
          SPGLib dataset - a dictionary, containing many informations about the symmetry and
          spacegroup of the structure. If it is not provided, and ``spacegroup`` is not ``None``
          the spacegroup is recomputed to obtain the dataset, if the dataset is requested.
        """
        atoms = sprkkr_atoms.SPRKKRAtoms.promote_ase_atoms(atoms)
        self.atoms = atoms
        self._dataset = dataset
        self.symmetry = symmetry
        self._block = None

    def to_dict(self):
        return {'dataset': self._dataset, 'symmetry': self.symmetry }

    """ ASE require this form of name"""
    todict = to_dict

    def copy_for(self, atoms):
        """ Return copy of the object for the given atoms """
        out = copy.copy(self)
        out.atoms = atoms
        out._dataset = None
        return out

    def check_spacegroup_kinds(self):
        """ Returns True, if the spacegroup_kinds array match the sites info. """
        a = self.atoms
        if 'spacegroup_kinds' in a.arrays:
            return False
        kinds = {}
        sites = {}
        occ = a.info.get('occupancy', {})
        for site, kind in zip(a.sites, a.array('spacegroup_kinds')):
            st = site.site_type
            if kinds.setdefault(kind, st) != st or \
               sites.setdefault(st, kind) != kind or \
               str(kind) not in occ:
                  return False
        return True

    def update_spacegroup_kinds(self, if_required=False):
        """ Update the occupancy info and spagroup_kinds array """
        if if_required and not self.check_spacegroup_kinds():
            return
        sgi = np.empty(len(self.atoms), dtype=int)
        stypes = {}
        occ = {}
        num = 0
        for i, site in enumerate(self.atoms.sites):
            st = site.site_type
            if not st in stypes:
                sgi[i] = stypes[st] = num
                occ[str(num)]= st.occupation.to_dict(True)
                num+=1
            else:
                sgi[i] = stypes[st]

        self.atoms.set_array('spacegroup_kinds', sgi)
        self.atoms.info['occupancy'] = occ

    def __repr__(self):
        return f"Spacegroup {str(self.spacegroup_number() or 'None')}"

    def __str__(self):
        return self.__repr__()

    def spacegroup_number(self) -> Optional[int]:
        """
        Returns
        -------
        spacegroup
          Spacegroup number or None, if there is no spacegroup.
        """
        dataset = self.dataset
        return dataset.number if dataset else None

    @property
    def dataset(self)->Optional[Dict]:
        """ Return SpgLib dataset containing informations about symmetry,
        spacegroup, equivalence of sites etc... """
        if self._dataset is None:
            self.recompute()
        return self._dataset

    @property
    def equivalent_sites(self)->np.ndarray:
        """ Return numpy array, that tags the sites by its equivalence
        classes """
        self.atoms.sites
        return self.atoms.get_array('spacegroup_kinds')

    def recompute(self,
                  symmetry=None,
                  consider_old=True,
                  precision=1e-5, angular_precision=1e-5,
                  *,
                  atomic_kinds=None,
                  copy=False,
                  init=False,
                  update_info=True
                  ):
       """ Init the sites array: the array of Sites objects,
       that holds additional SPRKKR properties and informations about site-equivalence """
       if symmetry is not None:
           self.symmetry = symmetry
       if self._block:
           if init:
               self.__block.init=True
           if copy:
               self._block.copy=True
           if update_info:
               self._block.update_info=True
           if atomic_kinds:
               raise NotImplementedError("Calling recompute during _block with atomic_kinds specified is not implemented")
           if not consider_old:
               self.consider_old=False
           self._block.angular_precision = min(angular_precision, self._block.angular_precision)
           self._block.precision = min(precision, self._block.precision)
           self._block.do = True
           return

       atoms = self.atoms
       if init or not atoms.are_sites_inited:
           init = True
       else:
           old_sites = atoms.sites
       creator = SiteType.creator(atoms)
       if not init:
           stype_copier = site_type_copier(atoms) if copy else used_site_type_copier(atoms)

           def create(i):
               osite = old_sites[i]
               return stype_copier(osite) if osite else creator(i)
       else:
           create = creator

       sites = np.zeros(len(atoms), dtype=object)
       indexes = np.arange(len(atoms))

       def set_site(sites, i, stype):
           if sites[i]:
               sites[i].unregister()
           sites[i] = Site(stype)

       def solve_region(region, slice):
            ssites = sites[slice]
            if sum(region.pbc) == 3 and self.symmetry or region is atoms:
                kinds = None if atomic_kinds is None else atomic_kinds[slice]
                dataset = spglib_dataset(region, kinds,
                                         consider_old = consider_old,
                                         precision=precision,
                                         angular_precision=angular_precision,
                                         add = ssites if atoms.regions else None)
            else:
                dataset = False

            to_global = indexes[slice]
            if self.symmetry and dataset:
                uniq, index, umap = np.unique(dataset.equivalent_atoms, return_index=True, return_inverse=True)
                stypes = np.empty(len(uniq), dtype=object)
                for i, site in enumerate(index):
                    stypes[i] = create(to_global[site])
                for i in range(len(umap)):
                    stype = stypes[umap[i]]
                    set_site(ssites, i, stype)

            else:
                for i in range(len(region)):
                    stype = create(to_global[i])
                    set_site(ssites, i, stype)

            if region is atoms:
                self._dataset = dataset

       for r in atoms.regions.values():
           solve_region(r, r.slice)
       solve_region(atoms, slice(None))
       atoms.set_sites(sites, True, update=update_info)
       return self

    @contextmanager
    def block_updating(self, always_recompute=False, **kwargs):
        if not self._block:
            block = self._block = SpacegroupInfoBlock()
        else:
            block = self._block.raise_counter()
        yield
        if always_recompute:
            self.recompute(**kwargs)
        self._block = self._block.lower_counter()
        block.recompute(self)


class SpacegroupInfoBlock:

    def __init__(self):
        self.init=False
        self.copy=False
        self.update_info=False
        self.consider_old=True
        self.precision = self.angular_precision=float('Inf')
        self.do = False
        self.counter = 1

    def raise_counter(self):
        self.counter += 1

    def lower_counter(self):
        self.counter -= 1
        return self if self.counter else None

    def recompute(self, info):
        if self.do:
            info.recompute(
                init=self.init,
                copy=self.copy,
                update_info=self.update_info,
                consider_old = self.consider_old,
                precision = self.precision,
                angular_precision = self.angular_precision,
            )
