""" The site class define the properties of an atom. """

from .radial_meshes import Mesh
from .reference_systems import ReferenceSystem
from .occupations import Occupation
import numpy as np

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
      self._occupation = Occupation.to_occupation(occupation)
      self._occupation._site = self
      self.atoms = atoms

  def copy(self):
      site = Site(self.atoms, self.occupation.copy(), self.reference_system.copy(), self.mesh.copy())
      site._occupation._site = site
      return site

  @property
  def occupation(self):
      """
      Creates the occupation, if not exists, according to the ASE atoms object
      """
      if self._occupation is None:
         ids = np.where(sites == self)[0]
         if not ids:
            raise ValueError('This atomic site is not from the provided Atoms object')
         an = atoms.get_atomic_numbers()
         oc = atoms.info.get('occupancy', {})
         for i in ids:
             if i in oc and oc[i]:
                self._occupation = Occupation(oc[i])
                return self._occupation
         for i in ids:
             if an[i]:
                self._occupation = Occupation(an[i])
                return self._occupation
         for i in ids:
             if atoms.symbols[i]:
                self._occupation = Occupation(atoms.symbols[i])
                return self._occupation
         raise ValueError('Unkwnown atom')

      return self._occupation

  def reset(self):
      self.mesh = Mesh.default()

  @occupation.setter
  def occupation(self, x):
      self._occupation = Occupation.to_occupation(x)
      self.update_atoms()

  @property
  def primary_symbol(self):
      return self.occupation.primary_symbol

  def update_atoms(self):
      index = self.atoms.sites == self
      an = self.atoms.get_atomic_numbers()
      an[ self.atoms.sites == self ] = self.occupation.primary_atomic_number
      self.atoms.set_atomic_numbers(an)
      ids = np.where(index)[0]
      occ = self.atoms.info.get('occupancy', {})
      for i in ids:
          occ[i] = self._occupation.as_dict
      self.atoms.info['occupancy'] = occ

  def __str__(self):
      return f"Site:{self.occupation}"

  def __repr__(self):
      return f"Site:{self.occupation}"
