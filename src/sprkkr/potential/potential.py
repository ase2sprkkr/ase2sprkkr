from ..common.conf_containers import RootConfContainer
from ..ase.atoms_io_data import AtomsIOData
from .definitions import potential

class Potential(RootConfContainer, AtomsWrapper):

  def __init__(self, atoms=None, definition=None):
      if definition is None:
         definition = definitions.potential
      super().__init__(definition=None)
      self._atoms_io_data = AtomsIOData.__init__(self, atoms)

  @property
  def atoms_io_data(self):
      return self._atoms_io_data

  @property
  def atoms(self):
      return self.atoms_io_data.atoms

  @atoms.setter
  def atoms(self, atoms):
      self.atoms_io_data.atoms = atoms

