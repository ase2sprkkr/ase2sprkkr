import cython
import ase
import numpy as np
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ase.spacegroup.spacegroup import Spacegroup
from ..spglib import SpacegroupInfo

cdef extern from "symmetry.h":
  cdef void find_symmetry_(int* n_operations, int* operations, int* ln, char* spacegroup, double* cell, double* angles, int* n, double* positions, cython.bint *align, double* magnetic, int* verbose)

_magnetic = np.zeros(3)

def find_symmetry_ex(spacegroup,
                  cython.double[::1] cell not None,
                  cython.double[::1] angles not None,
                  cython.int n,
                  cython.double[:,::1] positions not None,       #equivalent sites only!!!
                  cython.bint align=False,
                  cython.double[::1] magnetic_direction not None = _magnetic,
                  cython.bint verbose=False):
    """ Find point symmetry operations for a given cell,
    in the numbering used in xband """
    if isinstance(spacegroup, Spacegroup):
        spacegroup = spacegroup.no
    elif isinstance(spacegroup, SpacegroupInfo):
        spacegroup = spacegroup.number()
    spacegroup = str(spacegroup)
    cdef bytes sg = spacegroup.encode("utf8")
    cdef int ln = len(sg)
    cdef int _align = 1 if align else 0
    out = np.empty((2, 64), dtype=np.dtype('i'))
    cdef int[:,::1] _out = out
    cdef int n_out
    cdef int _verbose = 1 if verbose else 0
    find_symmetry_(&n_out, &_out[0,0], &ln, sg,
                   &cell[0], &angles[0], &n, &positions[0,0],
                   &align, &magnetic_direction[0], &_verbose)
    if n_out < 0:
        raise Exception("Failed to determine the symetry")
    out = out[:,:n_out]
    return out

def find_symmetry(atoms: ase.Atoms, align=False, verbose=False):
    SPRKKRAtoms.promote_ase_atoms(atoms)
    cp = atoms.cell.cellpar()
    es = atoms.spacegroup_info.equivalent_sites
    uniq = es.unique_indexes()

    """ Find point symmetry operations for a given cell,
    in the numbering used in xband """
    return find_symmetry_ex(
        atoms.spacegroup_info or 0,
        cp[:3],
        cp[3:],
        len(uniq),
        atoms.positions[uniq],
        align=align,
        verbose=verbose)
