from ase import Atoms
from ..common.unique_values import UniqueValuesMapping
import spglib
from ase.spacegroup import Spacegroup
import numpy as np
from ..sprkkr.sites import Site
from ..common.misc import numpy_index

class SprKkrAtoms(Atoms):
   """ ASE Atoms object extended by the data necessary for SPR-KKR calculations """

   @staticmethod
   def promote_ase_atoms(obj):
       """ Convert ASE Atoms object to the one usable by SprKkr.
           For the case of the usability it is a bit ugly hack: The __class__ attribute
           is replaced so the extra methods and properties of the objects will
           be available.
       """
       if obj and not isinstance(obj, SprKkrAtoms):
          if obj.__class__ is Atoms:
             obj.__class__ = SprKkrAtoms
          else:
             if not isinstance(obj, Atoms):
                raise(f'Can not promote class {obj} of class {obj.__class__} to {SprKkrAtoms}')

             class SprKKrAtomsEx(obj.__class__, SprKkrAtoms):
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
          raise RuntimeError('Spacegroup not found')
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
          mapping = UniqueValuesMapping(self.get_atomic_numbers())
          if consider_old and self._unique_sites:
            mapping = mapping.merge(self._unique_sites)
          if atomic_numbers:
            mapping = UniqueValuesMapping(atomic_numbers)

          spacegroup = self.compute_spacegroup_for_atomic_numbers(mapping.mapping, symprec=symprec)

        self.info['spacegroup'] = spacegroup

        tags = spacegroup.tag_sites(self.get_scaled_positions())
        sites = np.empty(len(tags), dtype=object)

        uniq, umap = np.unique(tags, return_inverse = True)
        for i in range(len(uniq)):
             index = umap == i
             if self._unique_sites is not None:
                #first non-none of the given index
                unique = next(filter(None, self._unique_sites[index]), None)
             else:
                unique = None
             sites[index] = unique or Site(self, self.symbols[ numpy_index(umap,i)])
        self.sites = sites

   @property
   def sites(self):
       if self._unique_sites is None:
          self.compute_sites_symmetry()
       return self._unique_sites

   @sites.setter
   def sites(self, v):
       an = np.zeros(len(v), dtype= int)
       for i,j in enumerate(v):
           an[i] = j.occupation.primary_atomic_number()
       self.set_atomic_numbers(an)
       self._unique_sites = v

   @property
   def potential(self):
       if self._potential is None:
          self._potential
