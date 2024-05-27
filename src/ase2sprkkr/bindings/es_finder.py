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
from ..physics.winger_seitz_radii import winger_seitz_radii
from typing import Dict, Union
from .empty_spheres import EmptySpheresResult


def empty_spheres(atoms: Atoms, *,
                  overlap_matrix:Union[float,np.ndarray]=0.18,
                  radii_ratios_map: Dict[str, float]=None,
                  max_es_overlap:float = 0.24,  # Maximum overlap of ES
                  adjust_overlap:float = 0.28,  # Overlap that will be adjusted to max. overlap
                  min_es_radius: float = 0.2,  # Min. accepted sphere radius
                  max_es_radius :float = 1.0,  # Max. accepted sphere radius
                  symmetrize_threshold: float = 0.7,  # Threshold for overlap when symmetrization is used
                  max_iterations: int = 100,  # Number of iterations for the sphere search
                  grid: np.ndarray = np.array([[48, 0, 0], [0, 48, 0], [0, 0, 48]]),
                  verbosity: int = 0,  # Verbosity of output
  ):
  """
  Compute the best coverage of the primitive cell with spheres.

  The function computes:
   - the best radii of atomic-sites spheres
   - the positions of 'vacuum pseudoatom' spheres and their radii to be added to the structure

  Parameters
  ----------
  atoms
    The structure, that will be filled in with empty spheres

  extend
    If True, the atoms will be extended by the empty spheres sites

  overlap_matrix
    #TODO

  radii_ratios_map
    Dict, that defines the ratios of radii of used chemical elements.
    It is dimensionless value (it is not the sizes, only
    the ratios of the given values, that matter).
    The default walues are :mod:`Winger-Seitz radii<ase2sprkkr.physics.winger_seitz_radii>`.
    If only some of the elements are provided, they are supplemented with the default values.

  max_es_overlap
    Max allowed overlap of empty spheres.

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

  params = dict(locals())
  for i in ('atoms', 'overlap_matrix', 'radii_ratios_map'):
      del params[i]

  if radii_ratios_map:
    radii_ratios_map = { **winger_seitz_radii, **radii_ratios_map }
  else:
    radii_ratios_map = winger_seitz_radii

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
  for s in map(str, pmg_structure.species):
    if not s in radii_ratios_map:
      radii_ratios_map[s] = 1.0
  structure = StructureAdapter(pmg_structure, radii_ratios_map, overlap_matrix)

  params = Parameters(**params)

  from ..common.no_output import NoOutput
  with NoOutput(suppress=verbosity<=0):
    res = run_finder(params, structure, symmetry)

  return EmptySpheresResult(res.es_positions, res.es_radii @ atoms.cell)
