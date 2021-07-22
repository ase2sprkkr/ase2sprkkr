""" This file contains SPRKKRAtoms - an enhanced version of Atoms to be used
with SPRKKR """


from ase import Atoms
from ..common.unique_values import UniqueValuesMapping
import spglib
from ase.spacegroup import Spacegroup
import numpy as np
from ..sprkkr.sites import Site
from ..common.misc import numpy_index

class SPRKKRAtoms(Atoms):
   """ ASE Atoms object extended by the data necessary for SPR-KKR calculations """

   @staticmethod
   def promote_ase_atoms(obj):
       """ Convert ASE Atoms object to the one usable by SPRKKR.
           For the case of the usability it is a bit ugly hack: The __class__ attribute
           is replaced so the extra methods and properties of the objects will
           be available.
       """
       if obj and not isinstance(obj, SPRKKRAtoms):
          if obj.__class__ is Atoms:
             obj.__class__ = SPRKKRAtoms
          else:
             if not isinstance(obj, Atoms):
                raise(f'Can not promote class {obj} of class {obj.__class__} to {SPRKKRAtoms}')

             class SprKKrAtomsEx(obj.__class__, SPRKKRAtoms):
                   pass
             obj.__class__ = SprKKrAtomsEx

          obj._init()
       return obj

   def __init__(self, *args, **kwargs):
       self._init()
       super().__init__(*args, **kwargs)

   def _init(self):
       """ The initialization of the additional (not-in-ASE) properties. To be used
       by constructor and by promote_ase_atoms"""
       self._unique_sites = None
       self._potential = None

   def compute_spacegroup_for_atomic_numbers(self, atomic_numbers=None, symprec=1e-5):
       """ Return spacegroup that suits to the atoms' cell structure and to the given
           atomic_numbers (not necessary the real ones, they can be just ''labels'').
       """

       atomic_numbers = atomic_numbers if atomic_numbers is not None else self.get_atomic_numbers()
       sg = spglib.get_spacegroup((self.get_cell(),
                             self.get_scaled_positions(),
                             atomic_numbers),
                             symprec=symprec)
       if sg is None:
           return None
       sg_no = int(sg[sg.find('(') + 1:sg.find(')')])
       spacegroup = Spacegroup(sg_no)
       return spacegroup

   def compute_sites_symmetry(self, spacegroup=None, atomic_numbers=None, consider_old=False, symprec=1e-5):
        """ SPRKKR has some properties shared by all by-symmetry-equal sites.
           This method initializes _sites property, that hold these properties:
           makes identical all atoms on the "symmetry identical positions with
           the same atomic number.

           It is called automatically when the sites property is firstly accessed.

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
              If True, and _unique_sites is not None, the non-symmetry-equivalent sites won't
              be equivalent in the newly computed symmetry.

            symprec: float
              A threshold for spatial error for symmetry computing. See spglib.get_spacegroup
        """
        if not spacegroup:
          if atomic_numbers:
            mapping = UniqueValuesMapping(atomic_numbers)
          else:
            mapping = UniqueValuesMapping(self.get_atomic_numbers())
            if consider_old and self._unique_sites:
              mapping = mapping.merge(self._unique_sites)
          spacegroup = self.compute_spacegroup_for_atomic_numbers(mapping.mapping, symprec=symprec)

        self.info['spacegroup'] = spacegroup

        occupation = self.info.get('occupancy', {})
        if spacegroup:
            tags = spacegroup.tag_sites(self.get_scaled_positions())
            sites = np.empty(len(tags), dtype=object)

            uniq, umap = np.unique(tags, return_inverse = True)
            used = set()
            for i in range(len(uniq)):
                 index = umap == i
                 if self._unique_sites is not None:
                    #first non-none of the given index
                    possible =  (i for i in self._unique_sites[index])
                    site = next(filter(None, possible), None)
                    if site in used:
                       site = site.copy()
                    else:
                       used.add(site)
                 else:
                    site = None
                 if not site:
                    symbol = self.symbols[ numpy_index(umap,i)]
                    for ai in np.where(index)[0]:
                        if ai in occupation and occupation[ai]:
                           symbol = occupation[ai]
                    site = Site(self, symbol)
                 sites[index] = site

        else:
            sites = np.empty(len(self), dtype=object)
            used = set()
            for i in range(len(self)):
                if self._unique_sites is not None:
                    site=self._unique_sites[i]
                    if site in used:
                       site = site.copy()
                    else:
                       used.add(site)
                else:
                    symbol = occupation[i] if i in occupation and occupation[i] else \
                             self.symbols[i]
                    site = Site(self, symbol)
                sites[i] = site

        self.sites = sites

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
       if self._unique_sites is None:
          self.compute_sites_symmetry()
       return self._unique_sites

   @sites.setter
   def sites(self, v):
       """ Set the sites property and update all other dependent
       properties (symbols, occupancy) according to the sites """
       an = np.zeros(len(v), dtype= int)
       occ = {}
       for i,j in enumerate(v):
           occ[i] = j.occupation.as_dict
           an[i] = j.occupation.primary_atomic_number
       self.set_atomic_numbers(an)
       self.info['occupancy'] = occ
       self._unique_sites = v

   @property
   def potential(self):
       if self._potential is None:
          self._potential = potentials.Potential.from_atoms(self)
       return self._potential

#at the last - to avoid circular imports
from ..potential import potentials
