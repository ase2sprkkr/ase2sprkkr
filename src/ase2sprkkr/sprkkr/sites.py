""" The site class define the properties of an atom. """

from .radial_meshes import Mesh
from .radial import RadialPotential, RadialCharge
from .reference_systems import ReferenceSystem
from .moments import Moments
import numpy as np
from ..common.decorators import cached_property
import copy
from typing import Optional, Union


class SiteType:
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
      self.atoms = atoms
      self.reference_system = reference_system or ReferenceSystem.default()
      self._mesh = mesh or Mesh.default()
      self._potential = None
      self._charge = None
      self._moments = None
      self._occupation = Occupation.to_occupation(occupation, None)
      self._occupation._site = self
      self.sites = set()

  def register(self, site):
      self.sites.add(site)

  def unregister(self, site):
      self.sites.remove(site)

  @property
  def has_symmetry(self):
      return len(self.sites) > 1

  def break_symmetry(self):
      if not self.sites:
          return
      itr = iter(self.sites)
      for i in itr:
          i.break_symmetry()

  def _clear_data(self):
      self._potential = None
      self._charge = None
      self._moments = None

  @property
  def mesh(self):
      return self._mesh

  @mesh.setter
  def mesh(self, mesh):
      self._mesh = mesh
      if self._potential:
          self._potential = self._potential.for_mesh(mesh)
      if self._charge:
          self._charge = self._charge.for_mesh(mesh)

  @property
  def potential(self):
      """ The radial potential data of the site """
      return self._potential

  @potential.setter
  def potential(self, value):
      if value is None:
          self._potential = None
          return
      if not isinstance(value, RadialPotential):
          value = RadialPotential(value, self._mesh, self.primary_atomic_number)
      self._potential = value
      self.mesh = value.mesh

  @property
  def charge(self):
      """ The radial charge data of the site """
      return self._charge

  @charge.setter
  def charge(self, value):
      if value is None:
          self._charge = None
          return
      if not isinstance(value, RadialCharge):
          value = RadialCharge(value, self._mesh)
      self._charge = value
      self.mesh = value.mesh

  @property
  def moments(self):
      """ The moments data of the site """
      return self._moments

  @moments.setter
  def moments(self, value):
      if not isinstance(value, Moments):
          value = Moments(*value)
      self._moments = value

  def copy(self, atoms=False):
      """ Create a copy of the site. """
      site_type = SiteType(self.atoms, self.occupation.copy(),
                           self.reference_system.copy(), self.mesh.copy())
      site_type._occupation._site_type = site_type
      if site_type.potential is not None:
          site_type.potential = site_type.potential.copy()
      if site_type.charge is not None:
          site_type.charge = site_type.charge.copy()
      if atoms is not False:
          site_type.atoms = atoms
      return site_type

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
         atoms = self.atoms
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

  @occupation.setter
  def occupation(self, x):
      self._occupation = Occupation.to_occupation(x, self)
      self.update_atoms()

  def reset(self):
      """
      Set the properties of the site to the default.

      Currently, it resets the mesh.
      """
      self.mesh = Mesh.default()
      self.potential = None
      self.charge = None

  @property
  def primary_symbol(self):
      """ Symbol of the most important (most probable) chemical element present on the site. """
      return self.occupation.primary_symbol

  @property
  def primary_atomic_number(self):
      """ Atomic symbol of the most important (most probable) chemical element present on the site. """
      return self.occupation.primary_atomic_number

  def index(self):
      """ Return the the sites-array (of the owning atoms object) index for this site. """
      enum = enumerate(i.site_type for i in self.atoms.sites)
      return [i[0] for i in enum if i[1] is self]

  def update_atoms(self):
      """ Update atomic numbers and occupation according to the sites data. """
      index = self.index()
      if not len(index):
         return
      an = self.atoms.get_atomic_numbers()
      pan = self.occupation.primary_atomic_number
      an[ index ] = pan
      if self._potential:
          self._potential.z = pan
      self.atoms.set_atomic_numbers(an)
      occ = self.atoms.info.get('occupancy', {})
      for i in index:
          occ[i] = self._occupation.as_dict
      self.atoms.info['occupancy'] = occ

  def __str__(self):
      return f"SiteType:{self.occupation}"

  def __repr__(self):
      return f"SiteType:{self.occupation}"

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


class Site:

  @staticmethod
  def create(atoms, occupation, reference_system=None, mesh=None):
      site_type = SiteType(atoms, occupation, reference_system, mesh)
      return Site(site_type)

  @staticmethod
  def copy_sites(sites, atoms=False):
      cache = {}

      def site(x):
          st = x._site_type
          n_st = cache.get(st, None)
          if not n_st:
             cache[st] = n_st = st.copy(atoms=False)
          return x.copy(site_type = n_st)

      try:
          return np.fromiter((site(i) for i in sites), dtype=object, count=len(sites))
      except ValueError:  # python 3.7 compatibility (or some old numpy version?)
          out = np.empty_like(sites)
          for i,s in enumerate(sites):
              out[i] = s
          return out

  def __init__(self, site_type):
      assert isinstance(site_type, SiteType)
      self._site_type = None
      self.site_type = site_type

  def __del__(self):
      if self._site_type:
          self._site_type.unregister(self)

  @property
  def site_type(self):
      return self._site_type

  @site_type.setter
  def site_type(self, site_type):
      if self._site_type:
          self._site_type.unregister(self)
      self._site_type = site_type
      if site_type:
          self._site_type.register(self)

  @property
  def has_symmetry(self):
      return self._site_type.has_symmetry

  def break_symmetry(self):
      if self.has_symmetry:
          self.site_type = self.site_type.copy()

  @property
  def mesh(self):
      return self._site_type.mesh

  @mesh.setter
  def mesh(self, mesh):
      self._site_type.mesh = mesh

  @property
  def potential(self):
      return self._site_type.potential

  @potential.setter
  def potential(self, potential):
      self._site_type.potential = potential

  @property
  def charge(self):
      return self._site_type.charge

  @charge.setter
  def charge(self, charge):
      self._site_type.charge = charge

  @property
  def moments(self):
      return self._site_type.moments

  @moments.setter
  def moments(self, moments):
      self._site_type.moments = moments

  @property
  def occupation(self):
      return self._site_type.occupation

  @occupation.setter
  def occupation(self, occupation):
      self._site_type.occupation = occupation

  @property
  def reference_system(self):
      return self._site_type.reference_system

  @reference_system.setter
  def reference_system(self, reference_system):
      self._site_type.reference_system = reference_system

  @property
  def atomic_types(self):
      return self._site_type.atomic_types

  @property
  def primary_symbol(self):
      """ Symbol of the most important (most probable) chemical element present on the site. """
      return self._site_type.primary_symbol

  @property
  def primary_atomic_number(self):
      """ Atomic symbol of the most important (most probable) chemical element present on the site. """
      return self._site_type.primary_atomic_number

  @property
  def is_vacuum(self):
      return self._site_type.is_vacuum

  def __str__(self):
      return f"Site:{self.occupation}"

  def __repr__(self):
      return f"Site:{self.occupation}"

  @property
  def atoms(self):
      return self._site_type.atoms

  def copy(self, site_type: Optional[Union[bool,SiteType]]=False):
      """
      Copy the site, optionally setting the new site_type (possibly to None) or
      copy the exiting one.
      """
      out = copy.copy(self)
      if site_type:
          if not isinstance(site_type, SiteType):
              site_type = out._site_type.copy()
      if site_type is not False:
          out._site_type = site_type
      if out._site_type:
          out._site_type.register(out)
      return out

  def _clear_data(self):
      return self._site_type._clear_data()

  def reset(self):
      return self._site_type.reset()


from .occupations import Occupation  # NOQA: E402
