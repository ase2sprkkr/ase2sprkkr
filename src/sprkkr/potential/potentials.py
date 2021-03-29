from ..common.conf_containers import RootConfContainer
from .io_data import ReadIoData, WriteIoData
from ..ase.sprkkr_atoms import SprKkrAtoms

class Potential(RootConfContainer):

  def __init__(self, atoms=None, definition=None,set_to_atoms=True):
      if definition is None:
         from .definitions.potential import potential_definition as definition
      self.atoms = atoms
      super().__init__(definition)
      if set_to_atoms and atoms:
         atoms.potential = self

  def read_from_file(self, file, atoms=None):
      super().read_from_file(file)
      if atoms is not False:
         self.atoms = self.update_atoms(atoms)

  def update_atoms(self, atoms=None):
      """ Update the ASE object from the values contained in the sections
      of the potential """
      atoms = SprKkrAtoms.promote_ase_atoms(atoms or self.atoms)
      iodata = ReadIoData()
      readed = set()

      def do(section):
          nonlocal atoms
          if section._definition.name in readed:
             return
          readed.add(section._definition.name)
          for d in section._depends_on():
             if d in self:
                do(self[d])
          atoms = section._update_atoms(atoms, iodata) or atoms

      for i in self:
          do(i)

      return atoms


  def save_to_file(self, file, atoms=None):
      if atoms is not False:
         self.set_from_atoms(atoms or self.atoms)
      super().save_to_file(file)

  def set_from_atoms(self, atoms = None):
      """ Set the sections of the potential accoring to an ASE atoms object. """
      atoms = SprKkrAtoms.promote_ase_atoms(atoms or self.atoms)
      iodata = WriteIoData(self.atoms)
      for i in self:
          i._set_from_atoms(atoms, iodata)

  @staticmethod
  def from_file(filename):
      pd = definition.potential_definition
      return pd.read_from_file(filename)

  @classmethod
  def from_atoms(cls, atoms, set_to_atoms=True):
      pd = definition.potential_definition
      return cls(atoms = atoms, definition = pd, set_to_atoms=set_to_atoms)

#on the last line to avoid the circular import issues
from .definitions import potential as definition
