import cython
import ase
import numpy as np
import warnings
cimport numpy as np

from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ase.spacegroup.spacegroup import Spacegroup
from ..spglib import SpacegroupInfo
from ...common.subprocess import in_subprocess
from ase.units import Bohr

cdef extern from "symmetry.h":
  cdef void find_symmetry_(int* n_operations, int* operations, int* ln, char* spacegroup, double* cell, double* angles, double* latvec,
                           int* n, double* cpositions, int* natom, int* types, double* positions, int* align, double* magnetic,
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

    find_symmetry_(&n_out, &_out[0,0], &ln, sg,
                   &cell[0], &angles[0], &latvec[0,0],
                   &n, &cpositions[0,0], &natom,
                   &types[0], &positions[0,0],
                   &_align, &magnetic_direction[0], &_verbose)
    if n_out < 0:
        raise Exception("Failed to determine the symetry")
    out = out[:,:n_out]
    return out

def find_symmetry(atoms: ase.Atoms, align=False, verbose=False, subprocess=True, use_spacegroup=True):
    """ Find point symmetry operations for a given cell,
    in the numbering used in xband """
    SPRKKRAtoms.promote_ase_atoms(atoms)
    es = atoms.spacegroup_info.equivalent_sites
    types = es.normalized(dtype=np.int32)[0]
    uniq = es.unique_indexes()

    if use_spacegroup or True:
       sg = atoms.spacegroup_info.number() or 0
    else:
       sg = 0
    if sg == 0:
       sg = atoms.spacegroup_info.number()
       msg = "Finding empty_spheres without providing " \
             "spacegroup is not currently supported."
       if sg == 0:
           warnings.warn(msg)
           return np.empty(0,2)
       warnings.warn(msg + " Fallbacking to use spacegroup")

    cdef double to_bohr = 1. / Bohr
    cdef double alat = atoms.cell.get_bravais_lattice().a * to_bohr
    cdef double ratio = to_bohr / alat

    cp = atoms.cell.cellpar()
    cp[:3] *= to_bohr
    spos = atoms.get_scaled_positions()
    cpos = spos[uniq]
    vectors = atoms.cell[:] * ratio
    pos = atoms.positions * ratio
    pos = np.ascontiguousarray(pos)
    cpos = np.ascontiguousarray(cpos)


    if subprocess:
        return in_subprocess('find_symmetry_ex', __name__,(
              sg,
              cp[:3],
              cp[3:],
              vectors,
              len(cpos),
              cpos,
              len(pos),
              types,
              pos,
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
            pos,
            align=align,
            verbose=verbose)
