""" This module contains routines for building systems, which can be used for SPR-KKR calculations
- e.g. system with vacuum pseudoatoms, or 2D semiinfinite systems
"""

from ..ase.build import aperiodic_times, stack as _stack
from .atoms_region import AtomsRegion
from .sprkkr_atoms import SPRKKRAtoms
import math
import numpy as np
from ase import Atoms
from ase.build import surface

from numbers import Real
from typing import Union, List, Tuple, Optional,Dict

def semiinfinite_system(atoms:Atoms, repeat:Union[Tuple[Real,Real],Real], atoms2:Atoms=None,
                        hkl:Optional[Tuple[Real]]=None, hkl2:Optional[Tuple[Real]]=None,
                        axis:int=2):
    """ Build a semiinfinite system from one or two 3D periodic systems.

    If two systems are given, they have to have identical the first two
    lattice vectors. If one system is given, the second system is created
    as a copy of the first system, but with vacuum (pseudo)atoms on its sites.

    Parameters
    ----------
    atoms
      The left bulk region of the result.

    repeat
      How many times should be the outer regions repeated in the central (non-bulk)
      region. If only one number is given, it is considered as the number of repeating
      of the left bulk region in the central region. The number of repeating of the right
      region then will be determined such that it will have the same integer part as the of
      the given number and the decimal part will be complement (to one).

    atoms2
      The right bulk region of the result. If it is None,
      it is created from the left one by replacing the atoms for vacuum pseudoatoms

    hkl
      If not None, rotate the left atoms according the given Miller coordinates, first.

    hkl2
      If not None, rotate the right atoms according the given Miller coordinates, first.
      If it is None, and the atoms2 are None too, the hkl argument is used for rotating
      the atoms2 object.

    axis
      Along which axis build the system.
    """

    if isinstance(repeat, (int, float)):
       repeat=(repeat, math.floor(repeat) + math.ceil(repeat) - repeat)

    if hkl is not None:
       atoms1 = surface(atoms, hkl, 1)
    else:
       atoms1 = atoms

    if atoms2 is None:
       atoms2 = vacuum_like(atoms1 if hkl2 is None else atoms)
    if hkl2 is not None:
       atoms = surface(atoms, hkl, 1)

    catoms2 = aperiodic_times(atoms2, repeat[1], axis=axis, direction=-1)
    catoms = aperiodic_times(atoms1, repeat[0], axis=axis)
    catoms = _stack([catoms, catoms2], axis=axis)

    out = stack( {'left': atoms1,
                  'central': catoms,
                  'right': atoms2}, axis=axis)
    return out

def stack(atomses:Dict[str,Atoms], axis:int, *args, inherit_cell=True, **kwargs):
    """ Stack the atoms along given axis to a one object, creating the regions in the atoms.

    The function accepts list of the names, the list of the Atoms objects, and all the other parameters
    of the function :func:`ase2sprkkr.ase.build.stack`.
    """
    cnt = None
    out = _stack(atomses.values(), axis=axis, *args, **kwargs)
    last = list(atomses.items())[-1]

    if inherit_cell is True:
       inherit_cell = np.ones(3, dtype=bool)
       inherit_cell[axis]=False

    for name, atoms in atomses.items():
        if atoms is last:
            upto = None
        else:
            upto = (cnt or 0) + len(atoms)
        AtomsRegion.from_atoms(atoms, name, slice(cnt, upto), inherit_cell=inherit_cell, atoms=out)
        cnt=upto
    return out


def vacuum_like(atoms):
    """
    Creates a copy of atoms, filled with vacuum pseudoatoms
    """
    out = atoms.copy()
    SPRKKRAtoms.promote_ase_atoms(out)
    for site in out.sites:
        site.occupation = { 'Vc' : 1.0 }
    return out
