from ..common.conf_containers import RootConfContainer
from ase import Atoms

class AtomsWrapper:

  def __init__(self):
      self._atoms = None

  @property
  def atoms(self):
      if self._atoms is None:
         self._atoms = Atoms()
      return self._atoms

  @atoms.setter
  def atoms(self, atoms):
      self._atoms = atoms

class Potential(RootConfContainer, AtomsWrapper):

  def __init__(self, definition):
      super().__init__(definition)
      AtomsWrapper.__init__(self)

