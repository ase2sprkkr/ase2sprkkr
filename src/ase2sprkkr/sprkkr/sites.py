""" The site class define the properties of an atom. """

from .radial_meshes import Mesh
from .reference_systems import ReferenceSystem
import numpy as np
from ..common.decorators import cached_property

class Site:
  """
  Definition of an atomic site.
  (By symmetry) equivalent sites should share the same definition.
  However, two same (by their atomic number) atoms in a spatially
  different positions (i.e. not by symmetry equal) should not share
  the same property.
  """

  def __init__(self, atoms, occupation, reference_system=None, mesh=None):
      """
      Parameters
      ----------
      occupation: dict
          { AtomicType: fraction }
          The type of the atom

          Occupation of the site. If None is given, the default occupation
          of { T : 1.0 } is used, where T is determined by ASE atomic number
          corresponding to the site.

      atoms: SPRKKRAtoms
          The atoms, into which the Site belongs

      reference_system: ReferenceSystem
          Default reference system is used if None is given

      mesh: Mesh
          Default ExponentialMesh is used if None is given
      """
      self.mesh = mesh or Mesh.default()
      self.reference_system = reference_system or ReferenceSystem.default()
      self._occupation = Occupation.to_occupation(occupation, None)
      self._occupation._site = self
      self.atoms = atoms

  def copy(self):
      """ Create a copy of the site. """
      site = Site(self.atoms, self.occupation.copy(), self.reference_system.copy(), self.mesh.copy())
      site._occupation._site = site
      return site

  @property
  def occupation(self):
      """
      The method returns the `:class:Occupation<ase2sprkkr.sprkkr.occupations.Occupation>` - that represent
      the occupation of the site by atoms. The occupation captures the probability of occurence of the
      given atoms on the site. There can be only partial occupancy (the probability, that there is
      an atom on the site is not 100%).

      The method creates the occupation, if it not exists, according to the ASE atoms object (using occpancy and symbols)
      """
      if self._occupation is None:
         ids = self.index()
         if not len(ids):
            raise ValueError('This atomic site is not from the provided Atoms object')
         an = atoms.get_atomic_numbers()
         oc = atoms.info.get('occupancy', {})
         for i in ids:
             if i in oc and oc[i]:
                self._occupation = Occupation(oc[i], self)
                return self._occupation
         for i in ids:
             if an[i]:
                self._occupation = Occupation(an[i], self)
                return self._occupation
         for i in ids:
             if atoms.symbols[i]:
                self._occupation = Occupation(atoms.symbols[i], self)
                return self._occupation
         raise ValueError('Unkwnown atom')

      return self._occupation

  def reset(self):
      """
      Set the properties of the site to the default.

      Currently, it resets the mesh.
      """

      self.mesh = Mesh.default()

  @occupation.setter
  def occupation(self, x):
      self._occupation = Occupation.to_occupation(x, self)
      self.update_atoms()

  @property
  def primary_symbol(self):
      """ Symbol of the most important (most probable) chemical element present on the site. """
      return self.occupation.primary_symbol

  def index(self):
      """ Return the the sites-array (of the owning atoms object) index for this site. """
      index = self.atoms.sites == self
      return np.where(index)[0]

  def update_atoms(self):
      """ Update atomic numbers and occupation according to the sites data. """
      index = self.index()
      if not len(index):
         return
      an = self.atoms.get_atomic_numbers()
      an[ index ] = self.occupation.primary_atomic_number
      self.atoms.set_atomic_numbers(an)
      occ = self.atoms.info.get('occupancy', {})
      for i in index:
          occ[i] = self._occupation.as_dict
      self.atoms.info['occupancy'] = occ

  def __str__(self):
      return f"Site:{self.occupation}"

  def __repr__(self):
      return f"Site:{self.occupation}"

  def is_vacuum(self) -> bool:
      """ Is the site vacuum pseudoatom? """
      return len(self._occupation) == 1 and next(iter(self.occupation)).is_vacuum()

  @cached_property
  def atomic_types(self):
      """
      This method provides the access to the atomic types of the current site.

      Returns
      -------
      atomic_types: AtomicTypesLookup
          Object, that can be indexed either by integer - position in the occupation (ordered) dictionary,
          or string - chemical symbol. It returns the atomic type, corresponding to the current position.
      """

      class AtomicTypesLookup:
          def __getitem__(_self, name):
              return self.occupation.atomic_type[name]
          def __setitem__(_self, name, value):
              return self.occupation.replace_type(name, value)

      return AtomicTypesLookup()

  @staticmethod
  def copy_sites(sites):
      cache = {}
      def site(x):
          if not x in cache:
             cache[x] = x.copy()
          return cache[x]
      return np.array([site(i) for i in sites], dtype=object)

from .occupations import Occupation
