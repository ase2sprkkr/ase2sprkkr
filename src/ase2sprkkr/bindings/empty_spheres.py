from collections import namedtuple
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms
import numpy as np


class EmptySpheresResult(namedtuple('EmptySpheresResult', 'positions radii')):

    def __len__(self):
        return len(self.radii)


def empty_spheres(atoms, method='auto', **kwargs):
    """ Returns centres of empty spheres to add """
    if method == 'auto':
        if es_finder.is_enabled:
            method = 'es_finder'
    if method == 'es_finder':
        return es_finder.empty_spheres(atoms, **kwargs)
    else:
        return spheres.empty_spheres(atoms, **kwargs)


def empty_spheres_atoms(atoms, round_zero=True, **kwargs):
    """
    Update the structure of the (SPRKKR) ASE atoms, adding the empty
    spheres and updating the shpheres radii of the atomic sites, according
    to an :func:`empty_spheres` result.
    """
    res = empty_spheres(atoms, **kwargs)
    num = len(res.radii)
    if num == 0:
       return None

    empty = SPRKKRAtoms(symbols='X' * len(res.radii),
                        positions=res.positions,
                        pbc = atoms.pbc,
                        cell = atoms.cell,
                        symmetry = False)
    empty.positions = empty.get_scaled_positions(True) @ empty.cell
    if round_zero:
        empty.positions[np.abs(empty.positions) < 1e-15] = 0
    return empty
    # for i,radius in zip(empty.sites, res.radii):
    # not TO DO: set the radius of the sphere, since SPRKKR make it itself
    #    next(iter(i.occupation)).radius=radius


def add_empty_spheres(atoms, *, copy=False, **kwargs):
    empty = empty_spheres_atoms(atoms, **kwargs)
    if empty:
        if copy:
            atoms = atoms.copy()
        atoms+= empty
    return atoms


from . import es_finder                   # NOQA
from .xband import spheres  # NOQA
