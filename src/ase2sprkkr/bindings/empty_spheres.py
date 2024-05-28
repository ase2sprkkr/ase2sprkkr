from collections import namedtuple
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms

EmptySpheresResult = namedtuple('EmptySpheresResult', 'positions radii')


def empty_spheres(atoms, method='auto', **kwargs):
    """ Returns centres of empty spheres to add """
    if method == 'auto':
        if es_finder.is_enabled:
            method = 'esfinder'
    if method == 'esfinder':
        return es_finder.empty_spheres(atoms, **kwargs)
    else:
        return spheres.empty_spheres(atoms, **kwargs)


def empty_spheres_atoms(atoms, **kwargs):
    """
    Update the structure of the (SPRKKR) ASE atoms, adding the empty
    spheres and updating the shpheres radii of the atomic sites, according
    to an :func:`empty_spheres` result.
    """
    res = empty_spheres(atoms, **kwargs)
    num = len(res.radii)
    if num == 0:
       return atoms

    empty = SPRKKRAtoms(symbols='X' * len(res.radii),
                        positions=res.positions,
                        pbc = atoms.pbc,
                        cell = atoms.cell,
                        symmetry = False)
    empty.positions = empty.get_scaled_positions(True) @ empty.cell
    return empty
    # for i,radius in zip(empty.sites, res.radii):
    # not TO DO: set the radius of the sphere, since SPRKKR make it itself
    #    next(iter(i.occupation)).radius=radius


def add_empty_spheres(atoms, **kwargs):
    empty = empty_spheres_atoms(atoms, **kwargs)
    atoms+= empty


from . import es_finder                   # NOQA
from .xband_empty_spheres import spheres  # NOQA
