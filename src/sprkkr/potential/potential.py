from ..common.conf_containers import RootConfContainer
from ase import Atoms
from ase.lattice import BravaisLattice
from ase.spacegroup import get_spacegroup

class AtomsWrapper:

  def __init__(self):
      self._atoms = None
      self.clear_cache()

  def clear_cache(self):
      self._cell_spacegroup = None
      self._bravais_lattice = None

  def set_atoms_cell(self, cell):
      self.clear_cache()
      if isinstance(cell, BravaisLattice):
         self._bravais_lattice = cell
         cell = cell.tocell()
      self.atoms.set_cell(cell)

  @property
  def atoms(self):
      if self._atoms is None:
         self._atoms = Atoms()
      return self._atoms

  @atoms.setter
  def atoms(self, atoms):
      self._atoms = atoms
      self.clear_cache()

  @property
  def bravais_lattice(self):
      if self._bravais_lattice is None:
          self._bravais_lattice = self.atoms.cell.get_bravais_lattice()
      return self._bravais_lattice

  @property
  def cell_spacegroup(self):
      if self._cell_spacegroup is None:
          self._cell_spacegroup = get_spacegroup(self.atoms)
      return self._cell_spacegroup

class Potential(RootConfContainer, AtomsWrapper):

  def __init__(self, definition):
      super().__init__(definition)
      AtomsWrapper.__init__(self)

