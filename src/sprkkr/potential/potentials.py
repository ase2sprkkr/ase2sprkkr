from ..common.conf_containers import RootConfContainer
from ..ase.atoms_io_data import AtomsIOData

class Potential(RootConfContainer):

  def __init__(self, atoms=None, definition=None):
      if definition is None:
         from .definitions.potential import potential_definition as definition
      super().__init__(definition=definition)
      self._atoms_io_data = AtomsIOData(atoms)

  @property
  def atoms(self):
      return self._atoms_io_data.atoms

  @atoms.setter
  def atoms(self, atoms):
      self._atoms_io_data.atoms = atoms
