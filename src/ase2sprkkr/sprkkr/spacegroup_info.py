"""
Wrapper for spglib for computing symmetry of a primitive cell
"""

from __future__ import annotations
from typing import Optional, Dict
import numpy as np
import copy
from contextlib import contextmanager

from . import sprkkr_atoms
from ..bindings.spglib import symmetry_dataset, symmetry_site_permutations
from .site_equivalence import site_equivalence_tags, refine_partial_site_tags
from .sites import Site, SiteType
from ..common.unique_values import UniqueValuesMapping


def site_type_copier(atoms):
    """Return a function, that always copy the site_type"""

    def copy(site):
        return site.site_type.copy(atoms=atoms)

    return copy


def used_site_type_copier(atoms):
    """Return a function, that copy the site_type, if it is requested more than once"""
    used = set()

    def copy(site):
        site_type = site.site_type
        if site_type in used or site_type.atoms is not atoms:
            site_type = site_type.copy(atoms=atoms)
        used.add(site_type)
        return site_type

    return copy


class BaseSpacegroupInfo:
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
    def dataset(self) -> Optional[Dict]:
        """Return SpgLib dataset containing informations about symmetry,
        spacegroup, equivalence of sites etc..."""
        if self._dataset is None:
            self.recompute()
        return self._dataset

    def kto_kyda_table(self, symprec=1e-5):
        """
        Build KTO_KYDA-style table, that maps the atoms before-after symmetry operations.

        Parameters
        ----------
        atoms: ase.Atoms
          The atoms object

        Returns
        -------
        table: (len(rotations), len(atoms) ) np.ndarray of
            If table[i,j] = k, then the i-th symmetry operation maps the j-th atom
            into k-th atom
        """

        atoms = self.for_atoms
        dataset = self.dataset
        if not dataset or dataset.rotations is None:
            return None

        tags = UniqueValuesMapping.from_values(
            [site.site_type for site in atoms.sites]
        ).mapping
        return symmetry_site_permutations(
            atoms.get_cell(), atoms.get_scaled_positions(), tags,
            dataset.rotations, dataset.translations, symprec=symprec,
        )


class RegionSpacegroupInfo(BaseSpacegroupInfo):
    def __init__(self, spacegroup_info, region):
        self.info = spacegroup_info
        self._dataset = None
        self.for_atoms = region

    def recompute(self):
        self.info.recompute()
        return self


class SpacegroupInfo(BaseSpacegroupInfo):
    """Class, that carry information about spacegroup and symmetry of a structure"""

    def __init__(self, atoms: sprkkr_atoms.SPRKKRAtoms, symmetry=True, dataset: Optional[Dict] = None):
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
        self._regions = {}

    @property
    def for_atoms(self):
        return self.atoms

    def to_dict(self):
        return {"dataset": self._dataset, "symmetry": self.symmetry}

    def info_for_region(self, region):
        name = region.name
        if name not in self._regions:
            self._regions[name] = RegionSpacegroupInfo(self, region)
        return self._regions[name]

    """ ASE requires this form of name"""
    todict = to_dict

    def copy_for(self, atoms):
        """Return copy of the object for the given atoms"""
        out = copy.copy(self)
        out.atoms = atoms
        out._dataset = None
        return out

    def check_spacegroup_kinds(self):
        """Returns True, if the spacegroup_kinds array match the sites info."""
        a = self.atoms
        if "spacegroup_kinds" not in a.arrays:
            return False
        kinds = {}
        sites = {}
        occ = a.info.get("occupancy", {})
        for site, kind in zip(a.sites, a.arrays["spacegroup_kinds"]):
            st = site.site_type
            key = str(int(kind))
            if (
                kinds.setdefault(kind, st) != st
                or sites.setdefault(st, kind) != kind
                or key not in occ
                or occ[key] != st.occupation.to_dict(True)
            ):
                return False
        return True

    def update_spacegroup_kinds(self, if_required=False, invalidate_spacegroup=False):
        """Update the occupancy info and spagroup_kinds array"""
        if invalidate_spacegroup:
            self._dataset = None
        if if_required and self.check_spacegroup_kinds():
            return
        sgi = np.empty(len(self.atoms), dtype=int)
        stypes = {}
        occ = {}
        num = 0
        for i, site in enumerate(self.atoms.sites):
            st = site.site_type
            if st not in stypes:
                sgi[i] = stypes[st] = num
                occ[str(num)] = st.occupation.to_dict(True)
                num += 1
            else:
                sgi[i] = stypes[st]

        self.atoms.set_array("spacegroup_kinds", sgi)
        self.atoms.info["occupancy"] = occ

    def __repr__(self):
        return f"Spacegroup {str(self.spacegroup_number() or 'None')}"

    def __str__(self):
        return self.__repr__()

    @property
    def equivalent_sites(self) -> np.ndarray:
        """Return numpy array, that tags the sites by its equivalence
        classes"""
        self.atoms.sites
        return self.atoms.get_array("spacegroup_kinds")

    def recompute(
        self,
        symmetry=None,
        consider_old=True,
        precision=1e-5,
        angular_precision=1e-5,
        *,
        atomic_kinds=None,
        copy=False,
        init=False,
        update_info=True,
    ):
        """Init the sites array: the array of Sites objects,
        that holds additional SPRKKR properties and informations about site-equivalence"""
        if symmetry is not None:
            self.symmetry = symmetry
        if self._block:
            if init:
                self._block.init = True
            if copy:
                self._block.copy = True
            if update_info:
                self._block.update_info = True
            if atomic_kinds is not None:
                raise NotImplementedError(
                    "Calling recompute during _block with atomic_kinds specified is not implemented"
                )
            if not consider_old:
                self._block.consider_old = False
            self._block.angular_precision = min(angular_precision, self._block.angular_precision)
            self._block.precision = min(precision, self._block.precision)
            self._block.do = True
            return

        atoms = self.atoms
        if init or not atoms.are_sites_inited():
            init = True
        else:
            # ASE set_array() overwrites an existing array in place. Keep the
            # previous Site instances for copying and unregistering below.
            old_sites = atoms.sites.copy()
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

        # Region membership is part of the site identity. Two sites from
        # different regions must never be merged by the global symmetry pass.
        memberships = [[] for _ in range(len(atoms))]
        for region_name, region in atoms.regions.items():
            for index in region.ids:
                memberships[int(index)].append(region_name)

        def set_site(target, i, stype, original=None):
            previous = target[i]
            site = original.copy(site_type=stype) if original else Site(stype)
            if previous:
                previous.unregister()
            target[i] = site

        def solve_region(region, slice):
            ssites = sites[slice]
            if sum(region.pbc) == 3 and self.symmetry or region is atoms:
                kinds = None if atomic_kinds is None else atomic_kinds[slice]
                tags = site_equivalence_tags(
                    region,
                    equivalence="site-types" if consider_old else "occupation",
                    atomic_kinds=kinds,
                    assigned_sites=ssites if atoms.regions else None,
                    region_memberships=(
                        memberships if region is atoms and atoms.regions else None
                    ),
                )

                def get_dataset(species):
                    try:
                        return symmetry_dataset(
                            region.get_cell(), region.get_scaled_positions(), species,
                            symprec=precision, angle_tolerance=angular_precision,
                        )
                    except ValueError:
                        return False

                dataset = get_dataset(tags)
                existing = (
                    region.arrays["sprkkr_sites_data"]
                    if "sprkkr_sites_data" in region.arrays else None
                )
                if (dataset is not False and consider_old and
                    existing is not None and any(not site for site in existing)):
                    refined = refine_partial_site_tags(tags, dataset, existing)
                    if not np.array_equal(tags, refined):
                        dataset = get_dataset(refined)
            else:
                dataset = False

            to_global = indexes[slice]

            def previous_site(local_index):
                if not init:
                    old = old_sites[to_global[local_index]]
                    if old:
                        return old
                current = ssites[local_index]
                return current if current else None

            def site_type_for(local_index):
                current = ssites[local_index]
                if region is atoms and atoms.regions and current:
                    # Preserve SiteTypes already established by a region.
                    return current.site_type
                return create(to_global[local_index])

            if self.symmetry and dataset:
                _, index, umap = np.unique(
                    dataset.equivalent_atoms, return_index=True, return_inverse=True
                )
                stypes = np.empty(len(index), dtype=object)
                representatives = [None] * len(index)
                for group, first in enumerate(index):
                    # Prefer a populated Site as the source of the shared
                    # SiteType, even when the first atom has no site yet.
                    source_index = first
                    if not init:
                        present = [
                            j for j in np.flatnonzero(umap == group)
                            if old_sites[to_global[j]]
                        ]
                        if present:
                            source_index = present[0]
                    representatives[group] = previous_site(source_index)
                    stypes[group] = site_type_for(source_index)
                for i, group in enumerate(umap):
                    original = previous_site(i) or representatives[group]
                    set_site(ssites, i, stypes[group], original)

            else:
                for i in range(len(region)):
                    set_site(ssites, i, site_type_for(i), previous_site(i))
            return dataset

        for r in atoms.regions.values():
            dataset = solve_region(r, r.slice)
            self.info_for_region(r)._dataset = dataset
        dataset = solve_region(atoms, slice(None))
        self._dataset = dataset
        atoms.set_sites(sites, True, update=update_info)
        if not init:
            # Rebuilt sites are copies, never the previous Site instances.
            # Unregister obsolete sites owned by these atoms. ASE slicing and
            # extending can temporarily borrow Sites from another Atoms object;
            # those must remain registered with their original owner. Repeated
            # arrays can also contain the same previous Site more than once.
            for previous in old_sites:
                if (previous and previous.site_type is not None
                        and previous.site_type.atoms is atoms):
                    previous.unregister()
        return self

    @contextmanager
    def block_updating(self, always_recompute=False, **kwargs):
        if self._block is None:
            block = self._block = SpacegroupInfoBlock()
        else:
            block = self._block
            block.raise_counter()
        try:
            yield
            if always_recompute:
                self.recompute(**kwargs)
        finally:
            self._block = block.lower_counter()
            if self._block is None:
                block.recompute(self)


class SpacegroupInfoBlock:
    """Class, for locking the update of spacegroup_info. If more operation are performed, the
    update is performed after the all operations are done."""

    def __init__(self):
        self.init = False
        self.copy = False
        self.update_info = False
        self.consider_old = True
        self.precision = self.angular_precision = float("Inf")
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
                consider_old=self.consider_old,
                precision=self.precision,
                angular_precision=self.angular_precision,
            )
