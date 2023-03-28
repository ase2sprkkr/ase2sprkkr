""" This module defines :class:`AtomsRegion`: a class,
that describe a symmetry of a region of Atoms object. """

from ase.cell import Cell
from ase import Atoms
import numpy as np
from typing import List, Union, Optional
from ..common.decorators import cached_property
from .sprkkr_atoms import SPRKKRAtoms

class AtomsRegion:
  """ AtomsRegion define a region of Atoms object,
  that can have its own cell and pbc.

  E.g. this system::

    |^^T^^^^^^^^^^^^^^^^^^^T^^|
    |  |                   |  |
    |__|___________________|__|

  has two semiinfinite 3D periodic regions and
  a central 2D-symmetric region. """

  def __init__(self, name:str, slice:slice,
               cell:Union[Cell, np.ndarray],
               pbc:List[Optional[bool]]=[True, True, True],
               inherit_cell:Union[bool,List[bool]]=False,
               atoms:Atoms=None):
      """
      Atoms region defines a spatial part of an Atoms object.
      The cell can has zero vectors - such vectors are taken from the parent
      atoms object.

      Parameters
      ----------
      slice
        Which sites (atoms) belong to the region

      cell
        The primitive cell

      pbc
        The peridodicity vector (see :class:`ase.Atoms`)

      inherit_cell
         Inherit the cell and the pcb from the parent atoms for given axes.
         Use of this argument is the same, as setting ``pcb`` to None and
         ``cell`` to ``[0,0,0]`` for the given axes.

      atoms
        The master object, of which the region is defined
      """
      self.name = name
      self.atoms = atoms
      self._slice = slice
      if inherit_cell is not False:
         cell = cell.copy()
         if inherit_cell is True:
             cell[:] = 0
             pbc = [ None, None, None ]
         else:
             pbc = np.array( pbc, dtype=object )
             cell[inherit_cell] = 0
             pbc[inherit_cell] = None
      self.incomplete_cell = Cell(cell)
      self.incomplete_pbc = pbc
      if atoms:
         self.set_atoms(atoms)

  @staticmethod
  def from_atoms(from_atoms: Atoms, name:str, slice:slice, inherit_cell=False, atoms:Atoms=None):
      """ Creates a region from the property of the given atoms object. The pbc and cell
      are taken from it. """
      return AtomsRegion(name, slice, from_atoms.cell, from_atoms.pbc, inherit_cell, atoms)

  @property
  def cell(self):
      """
      Return the cell.

      If the cell has some zero vectors, these are replaced with parent atoms cell vectors
      """
      cell = self.incomplete_cell.copy()
      if self.atoms is not None:
          for i,v in enumerate(cell):
              if not v.any():     #no nonzero
                  v[:] = self.atoms.cell[i]
      return cell

  @cell.setter
  def cell(self, cell):
      if not isinstance(cell, Cell):
          cell = np.asarray(cell)
          assert cell.shape == (3, 3)
      self.incomplete_cell = cell

  @property
  def pbc(self):
      """
      Return the pbc

      If the pbc has None items, these are replace with parent atoms pbc values
      """
      pbc = np.array(self.incomplete_pbc)
      if self.atoms is not None:
          for i,v in enumerate(pbc):
              if v is None:
                  pbc[i] = self.atoms.pbc[i]
      return pbc

  @pbc.setter
  def pbc(self, pbc):
      assert len(pbc) == 3
      self.incomplete_pbc = pbc

  @property
  def slice(self):
      return self._slice

  @slice.setter
  def slice(self, slice):
      self._slice = slice
      self._clear_cache

  def _clear_cache(self):
      for i in ('ids', 'set_of_ids'):
          if i in self.__dict__:
              delattr(self, i)

  def set_atoms(self, atoms, add=True):
      """ Set the master atoms - of which the region is described """
      if add:
         #atoms will call set_atoms again
         SPRKKRAtoms.promote_ase_atoms(atoms)
         atoms.add_region(self)
      else:
         self.atoms = atoms
         self._clear_cache()

  def create_atoms(self):
      """ Create an atoms object, that contains only the region """
      atoms = self.atoms[self.slice]
      atoms.cell = self.cell
      atoms.pbc = self.pbc

  @property
  def symbols(self):
      """ Return the symbols in the region """
      return self.atoms.symbols[self.slice]

  def __len__(self):
      return len(self.ids)

  def shared_ids_with(self, region):
      """ Return ids of the sites, that belongs to the both regions """
      if self.atoms != region.atoms:
         return false
      return self.set_of_ids.intersection(region.set_of_ids)

  @cached_property
  def ids(self):
      """ Returns ids of sites that belongs to the region """
      return np.arange(len(self.atoms))[self.slice]

  @cached_property
  def set_of_ids(self):
      """ Returns set of ids of sites that belongs to the region """
      return set(self.ids)

  @property
  def positions(self):
      return self.atoms.positions[self.slice]

  def only_vacuum_atoms(self):
      for site in self.atoms.sites[self.ids]:
          for at in site.occupation:
              if not at.is_vacuum():
                  return False
      return True
