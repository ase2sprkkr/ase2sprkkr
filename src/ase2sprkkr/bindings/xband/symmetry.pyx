import cython
import ase
import numpy as np
cimport numpy as np

from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ase.spacegroup.spacegroup import Spacegroup
from ..spglib import SpacegroupInfo
from ...common.subprocess import in_subprocess

cdef extern from "symmetry.h":
  cdef void find_symmetry_(int* n_operations, int* operations, int* ln, char* spacegroup, double* cell, double* angles, double* latvec,
                           int* n, double* cpositions, int* natom, int* types, double* positions, cython.bint *align, double* magnetic,
                           int* verbose)

_magnetic = np.zeros(3)

def find_symmetry_ex(spacegroup,
                  cython.double[::1] cell not None,
                  cython.double[::1] angles not None,
                  cython.double[:,::1] latvec not None,
                  cython.int n,
                  cython.double[:,::1] cpositions not None,       #equivalent sites only, scaled pos!!!
                  cython.int natom,
                  cython.int[::1] types not None,
                  cython.double[:,::1] positions not None,
                  cython.bint align=False,
                  cython.double[::1] magnetic_direction not None = _magnetic,
                  cython.bint verbose=False):
    """ Find point symmetry operations for a given cell,
    in the numbering used in xband """
    import traceback
    print()
    print()
    traceback.print_stack()

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
    cdef int i,j

    cdef np.ndarray latvec_ = np.empty((3,3))
    cdef cython.double[:,::1] _latvec_ = latvec_

    for i in range(3):
      for j in range(3):
        latvec_[i,j] = latvec[i,j] / cell[i]

    find_symmetry_(&n_out, &_out[0,0], &ln, sg,
                   &cell[0], &angles[0], &_latvec_[0,0],
                   &n, &cpositions[0,0], &natom,
                   &types[0], &positions[0,0],
                   &align, &magnetic_direction[0], &_verbose)
    if n_out < 0:
        raise Exception("Failed to determine the symetry")
    out = out[:,:n_out]
    return out

def find_symmetry(atoms: ase.Atoms, align=False, verbose=False, subprocess=True, use_spacegroup=True):
    """ Find point symmetry operations for a given cell,
    in the numbering used in xband """
    SPRKKRAtoms.promote_ase_atoms(atoms)
    cp = atoms.cell.cellpar()
    es = atoms.spacegroup_info.equivalent_sites
    types = es.normalized(dtype=np.int32)[0]
    uniq = es.unique_indexes()

    if use_spacegroup:
       sg = atoms.spacegroup_info.number() or 0
    else:
       sg = 0

    spos = atoms.get_scaled_positions()
    spos = np.ascontiguousarray(spos)
    cpos = spos[uniq]
    pos = atoms.positions
    vectors = atoms.cell[:]

    if subprocess:
        return in_subprocess('find_symmetry_ex', __name__,(
              sg,
              cp[:3],
              cp[3:],
              vectors,
              len(cpos),
              cpos,
              len(cpos),
              types,
              spos,
              ), {
                'align' :align,
                'verbose': verbose
              })
    else:
        return find_symmetry_ex(
            sg,
            cp[:3],
            cp[3:],
            vectors,
            len(cpos),
            cpos,
            len(pos),
            types,
            spos,
            align=align,
            verbose=verbose)
