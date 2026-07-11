from collections import namedtuple
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms
import numpy as np
import warnings


class EmptySpheresResult(namedtuple("EmptySpheresResult", "positions radii")):
    def __len__(self):
        return len(self.radii)


def empty_spheres(atoms, method="auto", **kwargs):
    """Return centres and radii of empty spheres to add.

    ``method='auto'`` uses the optional es_finder backend when both es_finder
    and Pymatgen are importable, otherwise it uses the ase2sprkkr "in-house"
    method taken from XBand.
    es_finder is not publicly available package. An explicitly requested but
    unavailable es_finder backend also falls back to the in-house implementation,
    after emitting a warning.
    """
    es_finder_available = es_finder.is_enabled

    if method == "auto":
        method = "es_finder" if es_finder_available else "inhouse"

    if method == "es_finder":
        if not es_finder_available:
            warnings.warn(
                "The es_finder method requires both es_finder and pymatgen, "
                "but at least one of them is not installed. Falling back to "
                "the inhouse empty-spheres method. Install optional support "
                'for pymatgen with: python -m pip install "ase2sprkkr[es_finder]"; '
                "es_finder itself must also be installed.",
                RuntimeWarning,
                stacklevel=2,
            )
            method = "inhouse"
        else:
            return es_finder.empty_spheres(atoms, **kwargs)

    if method in ("inhouse", "xband"):
        return spheres.empty_spheres(atoms, **kwargs)

    raise ValueError("Unknown empty-spheres method {!r}; use 'auto', 'es_finder', or 'inhouse'".format(method))


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

    empty = SPRKKRAtoms(
        symbols="X" * len(res.radii), positions=res.positions, pbc=atoms.pbc, cell=atoms.cell, symmetry=False
    )
    if round_zero:
        empty.positions[np.abs(empty.positions) < 1e-15] = 0
    empty.positions = empty.get_scaled_positions(True) @ empty.cell
    return empty
    # for i,radius in zip(empty.sites, res.radii):
    # not TO DO: set the radius of the sphere, since SPRKKR make it itself
    #    next(iter(i.occupation)).radius=radius


def add_empty_spheres(atoms, *, copy=False, **kwargs):
    empty = empty_spheres_atoms(atoms, **kwargs)
    if empty:
        if copy:
            atoms = atoms.copy()
        atoms += empty
    return atoms


from . import es_finder  # NOQA
from .xband import spheres  # NOQA
