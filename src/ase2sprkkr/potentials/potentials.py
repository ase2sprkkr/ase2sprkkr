""" Potential is object that holds data form SPR-KKR potential file """

from ..sprkkr.configuration import ConfigurationFile
from ..sprkkr.io_data import ReadIoData, WriteIoData
from ..common.decorators import class_property, cache
from io import StringIO


class Potential(ConfigurationFile):
  """ It holds data form SPR-KKR potential file

  It, in addition to being a containers for their sections, can read/write
  its properties from/to an ASE atoms object."""

  def __init__(self, atoms=None, definition=None):
      if definition is None:
         from .definitions.potential import potential_definition as definition
      self._atoms = atoms
      self._complete = False
      super().__init__(definition)

  def read_from_file(self, file, atoms=None, allow_dangerous=False):
      super().read_from_file(file, allow_dangerous=allow_dangerous)
      self.make_complete()
      if atoms is not False:
         self._atoms = self.update_atoms(atoms or self._atoms)

  def make_complete(self):
      """ Call this function, if you set manually all the properties necessary
      to create the atoms object """
      self._complete = True

  @property
  def atoms(self):
      if self._atoms is None and self._complete:
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

  def set_from_atoms(self, atoms=None, io_data:WriteIoData=None):
      """ Set the sections' values of the potential according to the given ASE atoms object.

      Parameters
      ----------
      atoms
        The atoms object, from which the data will be set. If it is None, the ``atoms`` property
        of the potential (``self.atoms``) is used.

      io_data
        The additional (in the time of the creation frozen) state of the atoms object,
        that contains e.g. numbering of the sites, atomic types etc.
        If is not set, it is created from the atoms.
      """
      super().set_from_atoms(atoms or self._atoms, io_data)
      self.make_complete()

  @class_property
  @cache
  def potential_definition(cls):
      # import here to avoid circular import issues
      from .definitions import potential as definition
      return definition.potential_definition

  @staticmethod
  def from_file(filename, atoms=None, allow_dangerous=False):
      """ Create a potential from a given potential file. """
      pd = Potential.potential_definition
      return pd.read_from_file(filename, atoms=atoms, allow_dangerous=allow_dangerous)

  @classmethod
  def from_string(cls, string:str, atoms=None, allow_dangerous=False):
      """ Create a potential from a string containing a content of a potential file. """
      return cls.from_file(StringIO(string), atoms, allow_dangerous)

  @classmethod
  def from_atoms(cls, atoms):
      """ Create a potential, that describes the given atoms object. """
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

  def __repr__(self):
      return "SPRKKR POTENTIAL"

  def __str__(self):
      return "SPRKKR POTENTIAL"


# At last - to avoid circular import problem
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms  # NOQA: E402
