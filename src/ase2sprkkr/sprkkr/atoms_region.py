""" This module defines :class:`AtomsRegion`: a class,
that describe a symmetry of a region of Atoms object. """

from ase.cell import Cell
from ase import Atoms
import numpy as np
from typing import List
from ..common.decorators import cached_property

class AtomsRegion:
  """ AtomsRegion define a region of Atoms object,
  that can have its own cell and pbc.

  E.g. this system

  |^^T^^^^^^^^^^^^^^^^^^^T^^|
  |  |                   |  |
  |__|___________________|__|

  has two semiinfinite 3D periodic regions and
  a central 2D-symmetric region. """

  def __init__(self, name:str, slice:slice, cell:Cell, pbc:List=[True, True, True], atoms:Atoms=None):
      """
      Parameters
      ----------
      slice
        Which sites (atoms) belong to the region

      cell
        The primitive cell

      pbc
        The peridodicity vector (see :class:`ase.Atoms`)

      atoms
        The master object, of which the region is defined
      """
      self.name = name
      self.atoms = atoms
      self.slice = slice
      self.cell = cell if isinstance(cell, Cell) else Cell(cell)
      self.pbc = np.asarray(pbc, dtype=bool)
      if atoms:
         atoms.add_region(self)

  def set_atoms(self, atoms):
      """ Set the master atoms - of which the region is described """
      self.atoms = atoms
      for i in ('ids', 'set_of_ids'):
          if i in self.__dict__:
              delattr(self, i)

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
