from ase.lattice import BravaisLattice
from ase.spacegroup import get_spacegroup
from .sprkkr_atoms import SprKkrAtoms

class AtomsIOData:
  """ Class used as interface between Atoms object and potential sections.
      It allows to store informations during storing (to avoid multiple
      evaluation of some quantities) and loading (to gather all necessary
      information to create the Atoms) of potential file.

      Do not use data provided by this object, if you do not know what you do,
      they may be obsolete.
  """

  def __init__(self, atoms=None):
      self.atoms = atoms

  def set_atoms(self):
      """ Set the Atom object to be computed with"""
      self._atoms = SprKkrAtoms.convert_from_ase_atoms(atoms)
      self.clear_cache()

  def create_atoms(self, **kwargs):
      """ Create new Atom object to be computed with """
      self._atoms = SprKkrAtoms(cell = self._cell, **kwargs)
      self.clear_cache()

  def set_new_positions(self, positions):
      """ Set new positions. Creates a new Atoms object, 
          if the numbers of atoms differ.
      """
      if self._atoms and len(self._atoms) == len(self.positions):
         self._atoms.positions = positions
      else:
         self.create_atoms(positions = positions)

  def clear_cache(self):
      """ Clear cached results (used during reading or writing of the POT file """
      self._cell_spacegroup = None
      self._bravais_lattice = None
      self._cell = None

  @property
  def cell(self):
      """ Return the Cell object, either from the current Atoms object or the one prepared
          for creating of the new Atoms class (during reading of POT file) """
      if self._atoms:
          return self._atoms.cell
      return self._cell

  @cell.setter
  def cell(self, cell):
      if isinstance(cell, BravaisLattice):
         self._bravais_lattice = cell
         cell = cell.tocell()
      if self._atoms:
         self._atoms.cell = cell
      self._cell = cell

  @property
  def atoms(self):
      return self._atoms

  @atoms.setter
  def atoms(self, atoms):
      self.set_atoms(atoms)

  @property
  def bravais_lattice(self):
      if self._bravais_lattice is None:
          self._bravais_lattice = self.cell.get_bravais_lattice()
      return self._bravais_lattice

  @property
  def cell_spacegroup(self):
      if self._cell_spacegroup is None:
          self._cell_spacegroup = get_spacegroup(self.atoms)
      return self._cell_spacegroup
