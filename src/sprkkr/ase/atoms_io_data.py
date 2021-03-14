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

  def set_atoms(self, atoms):
      """ Set the Atom object to be computed with"""
      self._atoms = SprKkrAtoms.convert_from_ase_atoms(atoms)
      self.clear_cache()

  def create_atoms(self, **kwargs):
      """ Create new Atom object to be computed with """
      self._atoms = SprKkrAtoms(cell = self._cell, **kwargs)
      self.clear_cache()

  def set_atoms_positions(self, positions):
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
  def bravais_cell(self):
      """ Return the Cell object, that will be used for generating the POT file.
      It is the cell given by bravais lattice of the current atoms. So, it can
      and can not be the same as the atoms cell """
      if not self._cell:
         cell = self._atoms.cell
         self.bravais_lattice = cell.get_bravais_lattice()
      return self._cell

  @property
  def atoms(self):
      return self._atoms

  @atoms.setter
  def atoms(self, atoms):
      self.set_atoms(atoms)

  @property
  def bravais_lattice(self):
      """ Bravais_lattice, beware, still in atomic units, so scaling on
          output is required """
      if self._bravais_lattice is None:
          self._bravais_lattice = self._atoms.cell.get_bravais_lattice()
      return self._bravais_lattice

  @bravais_lattice.setter
  def bravais_lattice(self, bl):
      self._bravais_lattice = bl
      self._cell = bl.tocell() if bl else None

  @property
  def cell_spacegroup(self):
      if self._cell_spacegroup is None:
          self._cell_spacegroup = get_spacegroup(self.atoms)
      return self._cell_spacegroup
