""" This file contains SPRKKRAtoms - an enhanced version of Atoms to be used
with SPRKKR """


from ase import Atoms
from ..bindings.spglib import SpacegroupInfo
import numpy as np
from ..sprkkr.sites import Site, SiteType
from ..common.misc import numpy_index
from typing import List, Union


class SPRKKRAtoms(Atoms):
   """ ASE Atoms object extended by the data necessary for SPR-KKR calculations """

   sites_array_name = 'sprkkr_sites_data'

   @staticmethod
   def promote_ase_atoms(obj, symmetry=None):
       """ Convert ASE Atoms object to the one usable by SPRKKR.
           For the case of the usability it is a bit ugly hack: The __class__ attribute
           is replaced so the extra methods and properties of the objects will
           be available.

           Parameters
           ----------
           obj: ase.Atoms
            The atoms object to be promoted to be used for SPRKKR calculations

           symmetry: bool or None
            The sites property of the resulting object will consider the symmetry of the structure.
            I.e., the by-symmetry-equal atomic sites will share the same sites object.
            Default None is the same as True, however it does not change the symmetry
            of the already promoted obj passed into the routine.
       """
       if obj and not isinstance(obj, SPRKKRAtoms):
          if obj.__class__ is Atoms:
             obj.__class__ = SPRKKRAtoms
          else:
             if not isinstance(obj, Atoms):
                raise ValueError(f'Can not promote class {obj} of class {obj.__class__} to {SPRKKRAtoms}')

             class SprKKrAtomsEx(obj.__class__, SPRKKRAtoms):
                   pass
             obj.__class__ = SprKKrAtomsEx

          obj._init(True if symmetry is None else symmetry)
          if obj.are_sites_inited():
             """ The sites array can be corrupted, since + and += operators of ase.Atoms
             do not handle it correctly. So: check it
             """
             for i in obj.sites:
                 if not i:
                      obj._init_sites(consider_old=True)
                      break

       else:
          if symmetry is not None:
             obj.symmetry = symmetry
       return obj

   def __init__(self, *args, symmetry=True, potential=None, **kwargs):
       """
       Creates SPRKKRAtoms atoms

       Parameters
       ----------
       *args: list
          The positionals arguments of ase.Atoms.__init__
       symmetry: bool
          The symmetry will be computed when the sites property will be initialized.
          I.e., the by-symmetry-equal atomic sites will share the same sites object.
       **kwargs: dict
          The named arguments of ase.Atoms.__init__
       """
       self._init(symmetry, potential)
       super().__init__(*args, **kwargs)

   def _init(self, symmetry=True, potential=None):
       """ The initialization of the additional (not-in-ASE) properties. To be used
       by constructor and by promote_ase_atoms"""
       self._symmetry = symmetry
       self._potential = potential
       self._regions = {}

   @property
   def regions(self):
       return self._regions

   def add_region(self, region):
       """
       Add a region of a given name
       """
       self._regions[region.name] = region
       region.set_atoms(self, False)

   def remove_region(self, name):
       del self._regions[name]

   def set_regions(self, regions:List):
       self._regions = { r.name: r for r in regions }
       for region in self._regions.values():
           region.set_atoms(self)

   @property
   def symmetry(self):
       """
       Whether the sites property is/will be generated using symmetry, i.e.
       whether the Sites objects in the sites property will be shared among
       symmetric atomic sites.
       """
       return self._symmetry

   @symmetry.setter
   def symmetry(self, value):
       """
       Recomputes the sites with enabled/disabled symmetry if the value of the property
       has changed.
       """
       if self._symmetry == value:
          return
       self._symmetry = value
       try:
         self.get_array(SPRKKRAtoms.sites_array_name)
         self._init_sites()
       except KeyError:
         pass

   def compute_sites_symmetry(self, consider_old=False, symmetry_precision=1e-5):
        """ SPRKKR has some properties shared by all by-symmetry-equal sites.
           This method initializes _sites property, that holds these properties:
           makes identical all the atoms on the "symmetry identical positions" with
           the same atomic number.

           The method is called automatically when the sites property is firstly accessed.
           The effect of the method is the nearly same as setting the symmetry property.
           However, setting the symmetry property on an 'already symmetrized' object has
           no effect, while this methods always recompute the sites property.

           Parameters
           ----------
           spacegroup: Spacegroup
              If not None, the given spacegroup is used for determining the symmetry,
              instead of the one determined by cell geometry.

           atomic_numbers: [ int ]
              Atomic numbers used to determine the spacegroup (if it is not given) to compute
              the symmetry. The atomic numbers can be ''virtual'', just to denote the equivalence
              of the sites.
              The array should have the same length as the number of atoms in the unit cell.
              If None, self.symbols are used.

            consider_old: bool
              If True, and self.sites is not None, the non-symmetry-equivalent sites won't
              be equivalent in the newly computed symmetry.

            symmetry_precision: float
              A threshold for spatial error for the symmetry computing. See spglib.get_spacegroup

        """
        self._symmetry = True
        return SPRKKRAtoms._init_sites(**locals())

   def _init_sites(self, consider_old=False, symmetry_precision=1e-5):
        """ See compute_sites_symmetry - this metod does just the same, but it does not set the symmetry property.

        All the hard work: finding the sites symmetry and distinguishing which sites are equivalent
        is done in SpacegroupInfo.from_atoms. Here are only the sites created according to the
        infomation obtained from the method.
        """
        sg_info = SpacegroupInfo.from_atoms(self, consider_old=consider_old, precision=symmetry_precision)

        sites = np.empty(len(self), dtype=object)
        occupation = self.info.get('occupancy', {})
        used = set()

        try:
           old_sites=self.get_array(SPRKKRAtoms.sites_array_name)
        except KeyError:
           old_sites=None

        def copy_site_type(site, i):
            site_type = site.site_type
            if site_type in used:
                site_type = site_type.copy(atoms=self)
            else:
                if site_type.atoms is not self:
                    site_type = site_type.copy(atoms=self)
                used.add(site_type)
            return site_type

        if sg_info.spacegroup:
          uniq, umap = np.unique(sg_info.equivalent_sites.mapping, return_inverse=True)
          stypes = np.empty(len(uniq), dtype=object)
          smap = np.empty(len(uniq), dtype=object)
          for i in uniq:
              smap[i] = set()
          for i,typ in enumerate(umap):
              smap[typ].add(i)

          stypes = np.empty(len(uniq), dtype=object)
          for i in uniq:
              if old_sites is not None:
                  # first non-none of the given index
                  possible = (j for j in smap[i] if old_sites[j])
                  try:
                     site = next(possible, None)
                     site_type = copy_site_type(old_sites[site], site) \
                                 if site is not None and old_sites[site] else None
                  except StopIteration:
                     site_type = None
              else:
                  site_type = None
              if not site_type:
                  symbol = self.symbols[ numpy_index(umap,i)]
                  occ = None
                  for ai in smap[i]:
                      occ = occupation.get(ai, None) or \
                            occupation.get(str(ai), None)
                      if occ:
                         symbol = occ
                         break
                  site_type = SiteType(self, symbol)
              stypes[i] = site_type
          for i in range(len(umap)):
              stype = stypes[umap[i]]
              if old_sites is not None and old_sites[i]:
                  if old_sites[i].atoms is self:
                      sites[i] = old_sites[i]
                      sites[i].site_type = stype
                  else:
                      sites[i] = old_sites[i].copy(site_type = stype)
              else:
                  sites[i] = Site(stype)
        else:
          for i in range(len(self)):
              if old_sites is not None and old_sites[i]:
                  site=old_sites[i]
                  stype = copy_site_type(site, i)
                  if site.atoms is self:
                     site.site_type = stype
                  else:
                     site = Site(stype)
              else:
                  symbol = occupation.get(i, None) or \
                           occupation.get(str(i), None) or \
                           self.symbols[i]
                  site = Site.create(self, symbol)
              sites[i] = site
        self.set_sites(sites, sg_info)
        return sites

   def cancel_sites_symmetry(self):
        """ Cancel the use of symmetry in the structure, i.e., makes the Site object
        uniqe (not shared) for each atomic site.

        Calling this method is nearly equivalent to setting the symmetry property
        to False, however, this method always recompute the sites object, while
        setting symmetry=False recomputes the sites property only if it was previously
        set to True.
        """
        self._symmetry = False
        self._init_sites()
        return self.sites

   @property
   def sites(self):
       """ The sites property holds all the information for the SPR-KKR package:
           atomic types (including number of semicore and valence electrons),
           occupancy, symmetries, meshes...
           Some of the properties are stored in the ASE atoms properties
           (e.g. occupancy, atomic symbol), however, ASE is not able to hold them
           all and/or to describe fully the SPR-KKR options; thus, these properties
           are hold in this array.

           The changes made on this array are reflected (as is possible)
           to the ASE properties, but the opposite does not hold - to reflect the changes
           in these properties please create a new Atoms object with given properties.
       """
       try:
          return self.get_array(SPRKKRAtoms.sites_array_name, copy=False)
       except KeyError:
          return self._init_sites()

   @sites.setter
   def sites(self, v):
       """ Set the sites property and updatei/clear all the other dependent
       properties (symbols, occupancy, spacegroup_info) according to the sites. """
       self.set_sites(v)

   def are_sites_inited(self):
       return SPRKKRAtoms.sites_array_name in self.arrays

   def set_sites(self, sites:np.ndarray, spacegroup_info:Union[SpacegroupInfo, bool, None]=None):
       """ Set the sites property and update all other dependent
       properties (symbols, occupancy) according to the sites.

       Unlike ``sites`` setter, this method  allow also set the spacegoup_info
       containing the computed informations about the symmetry.

       Parameters
       ----------
       sites
        The array of the :class:`ase2sprkkr.sprkkr.sites.Site` objects.

       sg_info
        Information, about the symmetry.
        If None is given, no info is available, the symmetry will be determined if needed.
        If False, there is no symmetry.
        If True, retain the current symmetry.
       """

       an = np.zeros(len(sites), dtype=int)
       occ = {}
       for i,j in enumerate(sites):
           occ[i] = j.occupation.as_dict
           an[i] = j.occupation.primary_atomic_number
       self.set_atomic_numbers(an)
       self.info['occupancy'] = occ

       if spacegroup_info is True:
          # retain the current symmetry
          pass
       elif spacegroup_info is False:
          self.info['spacegroup_info'] = SpacegroupInfo(None)
       elif spacegroup_info:
          self.info['spacegroup_info'] = spacegroup_info
       elif 'spacegroup_info' in self.info:
          # None have been given as spacegroup_info argument, delete any information about the symmetry
          del self.info['spacegroup_info']
       self.set_array(SPRKKRAtoms.sites_array_name, sites)

   @property
   def spacegroup_info(self):
       if not 'spacegroup_info' in self.info:
          self._init_sites()
       return self.info['spacegroup_info']

   def __getitem__(self, i):
       out = super().__getitem__(i)
       if isinstance(out, Atoms) and 'spacegroup_info' in out.info:
          del out.info['spacegroup_info']
       if self.are_sites_inited():
          out._init_sites(consider_old=True)
       return out

   @property
   def potential(self):
       if self._potential is None:
          self._potential = potentials.Potential.from_atoms(self)
       return self._potential

   def has_potential(self):
       return self._potential is not None

   def reset_sprkkr_potential(self):
       for i in self.sites:
           i.reset()
       if self._potential:
          self._potential.reset(update_atoms=False)
          self._potential.set_from_atoms()

   def extend(self, other):
       ln = len(self)
       out = super().extend(other)
       if self.symmetry:
          try:
              self.get_array(SPRKKRAtoms.sites_array_name, copy=False)
          except KeyError:
              return out
       if self.are_sites_inited():
          if self is other:
              print("COPY")
              for i in range(ln, ln + ln):
                  self.sites[i] = self.sites[i].copy()
          self._init_sites(consider_old=True)
       return out

   def copy(self):
       out = super().copy()
       SPRKKRAtoms.promote_ase_atoms(out)
       if self.are_sites_inited():
           out.set_sites(Site.copy_sites(out.sites, self), True)
       return out

   def __add__(self, other):
       if self.are_sites_inited():
           SPRKKRAtoms.promote_ase_atoms(other)
           other.sites
       out=super().__add__(other)
       return out

   def __iadd__(self, other):
       if self.are_sites_inited():
           SPRKKRAtoms.promote_ase_atoms(other)
           other.sites
       return super().__iadd__(other)


# at the last - to avoid circular imports
from ..potentials import potentials   # NOQA: E402
