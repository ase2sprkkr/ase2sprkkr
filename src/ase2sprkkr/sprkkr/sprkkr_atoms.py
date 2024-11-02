""" This file contains SPRKKRAtoms - an enhanced version of Atoms to be used
with SPRKKR """

from ase import Atoms
import numpy as np
from typing import List, Union
from .sites import Site
from .spacegroup_info import SpacegroupInfo


class SPRKKRAtoms(Atoms):
   """ ASE Atoms object extended by the data necessary for SPR-KKR calculations """

   sites_array_name = 'sprkkr_sites_data'

   @staticmethod
   def promote_ase_atoms(obj, symmetry=None, update_info=None):
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
           update_info:
            If True, always update spacegroup info
            If False, never
            If None, update it, if it seems to be inited, since it can be messed up by the
                     Atoms-manipulating routines
       """
       if obj is None:
          return None
       if isinstance(obj, SPRKKRAtoms):
          if symmetry is not None:
              obj.compute_sites_symmetry(symmetry)
       else:
          if symmetry is None:
              symmetry = True
          if obj.__class__ is Atoms:
             obj.__class__ = SPRKKRAtoms
          else:
             if not isinstance(obj, Atoms):
                raise ValueError(f'Can not promote class {obj} of class {obj.__class__} to {SPRKKRAtoms}')

             class SprKKrAtomsEx(obj.__class__, SPRKKRAtoms):
                   pass
             obj.__class__ = SprKKrAtomsEx

          obj._init(symmetry=symmetry)
          # if Atoms are summed with SPRKKRAtoms, the content of sites
          # can be wrong
          if update_info is not False:
              if update_info or obj.are_sites_inited():
                  obj.spacegroup_info.recompute()
       return obj

   def __init__(self, *args, potential=None, symmetry=True, **kwargs):
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
       super().__init__(*args, **kwargs)
       self._init(symmetry, potential)

   def _init(self, symmetry=True, potential=None):
       """ The initialization of the additional (not-in-ASE) properties. To be used
       by constructor and by promote_ase_atoms"""
       self._potential = potential
       self._regions = {}
       self.info['spacegroup_info'] = SpacegroupInfo(self, symmetry)

   @property
   def regions(self):
       return self._regions

   def add_region(self, region):
       """
       Add a region of a given name
       """
       self._regions[region.name] = region
       region.set_atoms(self, False)
       self.info['spacegroup_info'].recompute()

   def remove_region(self, name):
       del self._regions[name]

   def set_regions(self, regions:List):
       with self.info['spacegroup_info'].block_updating():
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
       return self.info['spacegroup_info'].symmetry

   @symmetry.setter
   def symmetry(self, value):
       """
       Recomputes the sites with enabled/disabled symmetry if the value of the property
       has changed.
       """
       if self.spacegroup_info.symmetry != value:
           self.compute_sites_symmetry(value)

   def compute_sites_symmetry(self, symmetry=None, consider_old=False, symmetry_precision=1e-5, angular_precision=1e-5):
       """ SPRKKR has some properties shared by all by-symmetry-equal sites.
       This method can be called to recompute the equivalent sites.
       """
       self.spacegroup_info.recompute(symmetry, consider_old, symmetry_precision, angular_precision)
       self.spacegroup_info.update_spacegroup_kinds()

   def breaks_sites_symmetry(self, *args, **kwargs):
        """ Cancel the use of symmetry in the structure, i.e., makes the Site object
        uniqe (not shared) for each atomic site.

        Calling this method is nearly equivalent to setting the symmetry property
        to False, however, this method always recompute the sites object, while
        setting symmetry=False recomputes the sites property only if it was previously
        set to True.
        """
        self.compute_sites_symmetry(False, *args, **kwargs)

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
       if not SPRKKRAtoms.sites_array_name in self.arrays:
          self.spacegroup_info.recompute(init=True)
       return self.arrays[SPRKKRAtoms.sites_array_name]

   @sites.setter
   def sites(self, v):
       """ Set the sites property and updatei/clear all the other dependent
       properties (symbols, occupancy, spacegroup_info) according to the sites. """
       self.set_sites(v)

   def are_sites_inited(self):
       return SPRKKRAtoms.sites_array_name in self.arrays

   def set_sites(self, sites:np.ndarray, spacegroup_info:Union[SpacegroupInfo, bool, None]=None, update:bool=True):
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

       update
        Pass False to prevent updating spacegroup_kinds and occupancy
       """
       an = np.zeros(len(sites), dtype=int)
       for i,j in enumerate(sites):
           an[i] = j.occupation.primary_atomic_number
       self.set_atomic_numbers(an)

       if spacegroup_info is True:
          # retain the current symmetry
          pass
       elif spacegroup_info is False:
          self.info['spacegroup_info'] = SpacegroupInfo(self, symmetry=False)
       elif spacegroup_info:
          self.info['spacegroup_info'] = spacegroup_info
       self.set_array(SPRKKRAtoms.sites_array_name, sites)
       if update:
          self.spacegroup_info.update_spacegroup_kinds()

   @property
   def spacegroup_info(self):
      return self.info['spacegroup_info']

   def __getitem__(self, i):
       out = super().__getitem__(i)
       if isinstance(out, Atoms) and 'spacegroup_info' in out.info:
          out.info['spacegroup_info'] = out.info['spacegroup_info'].copy_for(out)
       if self.are_sites_inited():
          out.info['spacegroup_info'].recompute(copy=True)
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

   def copy(self):
       out = super().copy()
       SPRKKRAtoms.promote_ase_atoms(out)
       if self.are_sites_inited():
           out.set_sites(Site.copy_sites(out.sites, out), self.spacegroup_info.copy_for(out))
       return out

   def _promote_if_other(self, other):
       if self.are_sites_inited():
           SPRKKRAtoms.promote_ase_atoms(other)
           other.sites
       elif isinstance(other, SPRKKRAtoms) and other.are_sites_inited():
           self.sites

   def __add__(self, other):
       return super().__add__(other)

   def __iadd__(self, other):
       return super().__iadd__(other)

   def extend(self, other):
       self._promote_if_other(other)
       ln = len(self)
       super().extend(other)
       if 'spacegroup_kinds' in self.arrays:
           sk = self.arrays['spacegroup_kinds']
           shift = np.max(sk[:ln])
           sk[ln:]+=shift
           if 'occupancy' in self.info or 'occupancy' in other.info:
               self.info.setdefault('occupancy', {})
               self.info['occupancy'].update({ str(int(i) + shift): v for i,v in other.info.get('occupancy', {}).items() })
       self.info['spacegroup_info']._dataset = None
       if self.are_sites_inited():
           sites = self.sites
           for i in range(ln, len(self)):
               sites[i] = sites[i].copy(self)


# at the last - to avoid circular imports
from ..potentials import potentials   # NOQA: E402
