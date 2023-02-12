from .sprkkr_atoms import SPRKKRAtoms
import math
from ase2sprkkr.ase.build import aperiodic_times, stack
from .atoms_region import AtomsRegion

def semiinfinite_system(atoms, repeat, atoms2=None,axis=2):
    """ Build a semiinfinite system from one or two 3D periodic systems.

    If two systems are given, they have to have identical the first two
    lattice vectors. If one system is given, the second system is created
    as a copy of the first system, but with vacuum (pseudo)atoms.
    """

    if isinstance(repeat, (int, float)):
       repeat=(repeat, math.floor(repeat) + math.ceil(repeat) - repeat)

    if atoms2 is None:
       atoms2 = vacuum_like(atoms)
    catoms2 = aperiodic_times(atoms2, repeat[1], axis=axis, direction=-1)
    catoms = aperiodic_times(atoms, repeat[0], axis=axis)
    catoms = stack([catoms, catoms2], axis=axis)
    out = stack([atoms, catoms, atoms2], axis=axis)
    la = len(atoms)
    la2 = len(atoms2)
    AtomsRegion('left', slice(0,la), atoms.cell, out)
    AtomsRegion('central', slice(la,-la2), atoms.cell, out)
    AtomsRegion('right', slice(-la2), atoms.cell, out)
    return out

def vacuum_like(atoms):
       out = SPRKKRAtoms(symbols='X' * len(atoms),
                            positions=atoms.positions,
                            cell=atoms.cell)
       for site in out.sites:
           site.occupation = { 'Vc' : 1.0 }
       return out
