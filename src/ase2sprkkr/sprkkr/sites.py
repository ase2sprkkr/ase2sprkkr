from .radial_meshes import Mesh
from .reference_systems import ReferenceSystem
from .occupations import Occupation

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
      Parameter
      ---------
      occupation: dict
          { AtomicType: fraction }
          The type of the atom

          Occupation of the site. If None is given, the default occupation
          of { T : 1.0 } is used, where T is determined by ASE atomic number
          corresponding to the site.

      atoms: SprKkrAtoms
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

  @property
  def occupation(self):
      """
      Makes sure, that occupation of the
      """
      if self._occupation is None:
         ids = self.atoms.sites.idsof(self)
         if not ids:
            raise ValueError('This atomic site is not from the provided Atoms object')
         an = atoms.get_atomic_numbers()
         for i in ids:
             if an[i]:
                self.occupation = Occupation(an[i])
                return self.occupation
         for i in ids:
             if atoms.symbols[i]:
                self.occupation = Occupation(atoms.symbols[i])
                return self.occupation
         raise ValueError('Unkwnown atom')

      return self._occupation

  def update_atoms(self):
      an = self.atoms.get_atomic_numbers()
      an[ self.atoms.sites == self ] = self.occupation.primary_atomic_number()
      self.atoms.set_atomic_numbers(an)

  def __str__(self):
      return f"Site:{self.occupation}"

  def __repr__(self):
      return f"Site:{self.occupation}"
