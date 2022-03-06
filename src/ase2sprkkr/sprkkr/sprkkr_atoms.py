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
   def promote_ase_atoms(obj, symmetry=None):
       """ Convert ASE Atoms object to the one usable by SPRKKR.
           For the case of the usability it is a bit ugly hack: The __class__ attribute
           is replaced so the extra methods and properties of the objects will
           be available.

           Parameters
           ----------
           obj: ase.Atoms
            The atoms object to be promoted to be used for SPRKKR calculations

           symmetry: boolean or None
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
                raise(f'Can not promote class {obj} of class {obj.__class__} to {SPRKKRAtoms}')

             class SprKKrAtomsEx(obj.__class__, SPRKKRAtoms):
                   pass
             obj.__class__ = SprKKrAtomsEx

          obj._init(True if symmetry is None else symmetry)
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
       symmetry: boolean
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
       self._unique_sites = None
       self._potential = potential
       self._symmetry = symmetry

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
       if self._unique_sites is not None:
          if value:
             self._compute_sites_symmetry()
          else:
             self._cancel_sites_symmetry()


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
           This method initializes _sites property, that holds these properties:
           makes identical all atoms on the "symmetry identical positions with
           the same atomic number.

           The method is called automatically when the sites property is firstly accessed.
           The effect of the method is the nearly same as setting the symmetry property.
           However, setting the symmetry property on an `already symmetrized' object has
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
              If True, and _unique_sites is not None, the non-symmetry-equivalent sites won't
              be equivalent in the newly computed symmetry.

            symprec: float
              A threshold for spatial error for the symmetry computing. See spglib.get_spacegroup

        """
        self._symmetry = True
        SPRKKRAtoms._compute_sites_symmetry(**locals())

   def _compute_sites_symmetry(self, spacegroup=None, atomic_numbers=None, consider_old=False, symprec=1e-5):
        """ See compute_sites_symmetry - this metod does just the same, but it does not set the symmetry property."""

        occupation = self.info.get('occupancy', {})
        if not spacegroup and self._symmetry:
          if atomic_numbers:
            mapping = UniqueValuesMapping(atomic_numbers)
          else:
            mapping = UniqueValuesMapping(self.get_atomic_numbers())
            if consider_old and self._unique_sites:
              mapping = mapping.merge(self._unique_sites)
          if occupation:
            def gen_occ():
                for i in range(len(mapping)):
                    val = occupation.get(i, None)
                    if val is None:
                        yield val
                    else:
                        yield tuple((k, val[k]) for k in val)
            mapping = mapping.merge(gen_occ())

          spacegroup = self.compute_spacegroup_for_atomic_numbers(mapping.mapping, symprec=symprec)

        self.info['spacegroup'] = spacegroup

        if not spacegroup:
            return self.cancel_sites_symmetry()

        tags = spacegroup.tag_sites(self.get_scaled_positions())
        mapping = mapping.merge( tags )
        tags = mapping.mapping

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
        self.sites = sites

   def cancel_sites_symmetry(self):
        """ Cancel the use of symmetry in the structure, i.e., makes the Site object
        uniqe (not shared) for each atomic site.

        Calling this method is nearly equivalent to the setting the symmetry property
        to False, however, this method always recompute the sites object, while
        setting symmetry=False recomputes the sites property only if it was previously
        set to False.
        """
        self._symmetry = False
        self._cancel_sites_symmetry()

   def _cancel_sites_symmetry(self):
        """ See cancel_sites_symmetry - this metod does just the same, but it does not set the symmetry property."""
        sites = np.empty(len(self), dtype=object)
        used = set()
        occupation = self.info.get('occupancy', {})
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
          self._compute_sites_symmetry()
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

   @potential.setter
   def potential(self, potential):
       self._potential = potential

   def reset_sprkkr_potential(self):
       for i in self.sites:
           i.reset()
       if self._potential:
          self._potential.reset(update_atoms = False)
          self._potential.set_from_atoms()



#at the last - to avoid circular imports
from ..potentials import potentials
