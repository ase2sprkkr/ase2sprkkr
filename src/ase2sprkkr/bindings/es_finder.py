"""
Binding for the es_finder package, that can determine the ideal positions
of empty spheres to fill the gaps in the primitive cell.
"""
from __future__ import annotations

import numpy as np

try:
  from es_finder.core.driver import run_finder
  from es_finder.core.parameters import Parameters
  from es_finder.core.symmetry import Symmetry
  from es_finder.adapter.pymatgen import StructureAdapter
  from pymatgen.io.ase import AseAtomsAdaptor
  import_error = None
  is_enabled = True
except ImportError as error:
  import_error = error
  run_finder = None
  is_enabled = False

from ase import Atoms
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms
from typing import Dict


def empty_spheres(atoms: Atoms,
                  extend: bool = False,

                  overlap_matrix:float|np.ndarray=0.14,
                  radii_ratios_map: Dict[str, float]=None,
                  max_es_overlap:float = 0.26,  # Maximum overlap of ES
                  adjust_overlap:float = 0.36,  # Overlap that will be adjusted to max. overlap
                  max_es_radius :float = 3.0,  # Max. accepted sphere radius
                  min_es_radius: float = 0.2,  # Min. accepted sphere radius
                  symmetrize_threshold: float = 0.8,  # Threshold for overlap when symmetrization is used
                  max_iterations: int = 2,  # Number of iterations for the sphere search
                  grid: np.ndarray = np.array([[4, 0, 0], [0, 192, 0], [0, 0, 192]]),
                  verbosity: int = 0,  # Verbosity of output
    ):
  """
  Compute "empty spheres" that optimally fills the empty space in the primitive cell

  Parameters
  ----------
  atoms
    The structure, that will be filled in with empty spheres

  extend
    If True, the atoms will be extended by the empty spheres sites

  overlap_matrix
    #TODO

  radii_ratios_map
    Dict, that defines the radii of given chemical elements
    #TODO - generate it automatically

  max_es_overlap
    Max allowed overlap of empty spheres

  adjust_overlap
    #TODO

  max_es_radius
    Maximum radius of empty shperes

  min_es_radius
    Minimum radius of empty shperes

  symetrize_threshold
    #TODO
    Threshold for overlap when symmetrization is used

  max_iterations
    Number of iterations for sphere finding algorithm

  grid
    #TODO

  verbosity
    Verbosity of the output of the algorithm

  """
  if run_finder is None:
     raise ImportError('Cannot import es_finder or pymatgen. Please install it, or set empty_spheres parameter to False') \
         from import_error

  params = locals()
  for i in ('atoms', 'overlap_matrix', 'radii_ratios_map', 'extend'):
      del params[i]

  if radii_ratios_map is None:
    radii_ratios_map = {}
  for s in atoms.symbols:
    if not s in radii_ratios_map:
      radii_ratios_map[s] = 1.0

  atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
  sym_dataset = atoms.spacegroup_info.dataset
  if sym_dataset is None:
     sym_dataset = {
         'translations' : np.empty((0,3)),
         'rotations' : np.empty((0,3,3)),
         'number' : 0,
         'hall_number' : 0,
         'hall' : '',
     }
  symmetry=Symmetry(
      sym_dataset['rotations'],
      sym_dataset['translations'],
      sym_dataset['number'],
      sym_dataset['hall_number'],
      sym_dataset['hall']
  )

  pmg_structure = AseAtomsAdaptor.get_structure(atoms)
  structure = StructureAdapter(pmg_structure, radii_ratios_map, overlap_matrix)

  params = Parameters(**params)

  res = run_finder(params, structure, symmetry)
  num = len(res.es_radii)
  if num == 0:
     return None

  empty = SPRKKRAtoms(symbols='X'*len(res.es_radii), symmetry = False)
  empty.set_positions(res.es_positions)
  for i,radius in zip(empty.sites, res.es_radii):
      next(iter(i.occupation)).radius=radius
  if extend:
     atoms.extend(empty)
  return empty
