""" Potential is object that holds data form SPR-KKR potential file """

from ..common.configuration_containers import RootConfigurationContainer
from .io_data import ReadIoData, WriteIoData
from ..common.misc import class_property, cache

class Potential(RootConfigurationContainer):
  """ It holds data form SPR-KKR potential file

  It, in addition to being a containers for their sections, can read/write
  its properties from/to an ASE atoms object."""

  def __init__(self, atoms=None, definition=None):
      if definition is None:
         from .definitions.potential import potential_definition as definition
      self._atoms = atoms
      self._complete = False
      super().__init__(definition)

  def read_from_file(self, file, atoms=None):
      super().read_from_file(file)
      self.make_complete()
      if atoms is not False:
         self._atoms = self.update_atoms(atoms or self._atoms)

  def make_complete(self):
      """ Call this function, if you set manually all the properties necessary
      to create the atoms object """
      self._complete = True

  @property
  def atoms(self):
      if not self._atoms and self._complete:
         self._atoms = self.update_atoms()
      return self._atoms

  @atoms.setter
  def atoms(self, atoms):
      self._atoms = atoms

  def update_atoms(self, atoms=None):
      """ Update the ASE object from the values contained in the sections
      of the potential """
      atoms = SPRKKRAtoms.promote_ase_atoms(atoms or self._atoms)
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

  def save_to_file(self, file, atoms=None, *, validate=True):
      if atoms is not False:
         self.set_from_atoms(atoms or self._atoms)
      self.make_complete()
      super().save_to_file(file, validate=validate)

  def set_from_atoms(self, atoms = None):
      """ Set the sections of the potential accoring to an ASE atoms object. """
      atoms = SPRKKRAtoms.promote_ase_atoms(atoms or self._atoms)
      iodata = WriteIoData(self.atoms)
      for i in self:
          i._set_from_atoms(atoms, iodata)

  @class_property
  @cache
  def potential_definition(cls):
      #import he to avoid the circular import issues
      from .definitions import potential as definition
      return definition.potential_definition

  @staticmethod
  def from_file(filename, atoms=None):
      pd = Potential.potential_definition
      return pd.read_from_file(filename, atoms=atoms)

  @classmethod
  def from_atoms(cls, atoms):
      pd = Potential.potential_definition
      return cls(atoms = atoms, definition = pd)

  def reset(self, update_atoms=True):
      copy = [i for i in self]
      for i in copy:
          if not i._definition.mandatory:
             i.reset()
      self.SCF_INFO.SCFSTATUS = 'START'
      if update_atoms:
          self.update_atoms()

#At last - to avoid circular import problem
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms
