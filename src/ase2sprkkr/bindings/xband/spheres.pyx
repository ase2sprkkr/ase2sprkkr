import cython
import numpy as np
cimport numpy as np
from ase import Atoms
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ...common.unique_values import UniqueValuesMapping
from .symmetry import find_symmetry
from ..empty_spheres import EmptySpheresResult
from ase.units import Bohr


cdef extern from "spheres.h":
  cdef void find_empty_spheres_(
              int* n,
              double* centres,
              double* radii,
              double* min_radius,
              double* max_radius,
              double* alat,
              double* cell,
              int* n,
              double* positions,
              int* atom_eq_class,
              int* n_classes,
              int* n_types,
              char* symbols,
              double* atomic_numbers,
              double* occupations,
              int* type_eq_class,
              int* n_symop,
              int* symop_number,
              int* symop_data,
              int* mesh,
              int* n_symmetry_operations,
              double *rotations,
              double *translations,
              int* verbose,
  );


def empty_spheres(
    atoms,
    double min_radius=0.65,
    double max_radius=2.,
    int[:,:] point_symmetry=None,
    verbose=False,
    int max_spheres=256,
    return_atom_rws=False,
    use_spacegroup=True,
    mesh=24,
    double[:,:,::1] rotations=None,
    double[:,::1]   translations=None
    ):

    if point_symmetry is None:
        point_symmetry = find_symmetry(atoms, use_spacegroup=use_spacegroup)
    else:
        SPRKKRAtoms.promote_ase_atoms(atoms)

    cdef double to_bohr = 1 / Bohr
    min_radius *= to_bohr
    max_radius *= to_bohr
    cdef double alat = atoms.cell.get_bravais_lattice().a * to_bohr
    cdef int n = len(atoms)

    es = atoms.spacegroup_info.equivalent_sites
    es = UniqueValuesMapping(es)
    ui = es.unique_indexes()
    cdef int[:] mapping = es.normalized(dtype=np.int32)[0]

    cdef int i
    cdef int n_types = sum(len(atoms.sites[i].occupation) for i in ui)
    cdef int n_classes = len(ui)
    cdef np.ndarray symbols = np.zeros((n_types, 4), dtype=np.byte)
    cdef char[:,:] _symbols = symbols
    cdef double[:] occupations = np.empty(n_types, np.double)
    cdef double[:] atomic_numbers = np.empty(n_types, np.double)
    cdef int[:] eq_classes = np.empty(n_types, np.int32)
    cdef int[3] _mesh = np.empty(3, np.int32)
    cdef int type_no = 0
    if isinstance(mesh, int):
      _mesh[0] = _mesh[1] = _mesh[2] = mesh
    else:
      _mesh[:] = mesh

    for i, index in enumerate(ui):
        site = atoms.sites[index]
        for typ, occ  in site.occupation.items():
            s = typ.symbol.encode('utf8')[:4]
            symbols[type_no, :len(s)] = memoryview(s)
            occupations[type_no] = occ
            atomic_numbers[type_no] = float(typ.atomic_number)
            eq_classes[type_no] = index + 1
            type_no += 1

    cdef int n_symops=point_symmetry.shape[1]
    cdef double ratio = to_bohr / alat
    cdef double[:,:] cell = atoms.cell[:] * ratio
    cdef double[:,:] positions = atoms.positions * ratio
    cdef int _verbose

    if isinstance(verbose, bool):
        _verbose = 1 if verbose else 0
    else:
        _verbose = verbose

    centres = np.empty((max_spheres, 3), dtype=np.double)
    radii = np.empty(max_spheres+n, dtype=np.double)

    cdef double[:,:] _centres = centres
    cdef double[:] _radii = radii

    cdef int n_symmetry_ops= -1 if (translations is None) else len(translations)
    assert n_symmetry_ops == -1 if (rotations is None) else len(rotations)
    if n_symmetry_ops > 0:
        assert translations.shape[1] == 3
        assert rotations.shape[1] == 3
        assert rotations.shape[2] == 3

    find_empty_spheres_(
                   &max_spheres,
                   &_centres[0,0],
                   &_radii[0],
                   &min_radius,
                   &max_radius,
                   &alat,
                   &cell[0,0],
                   &n,
                   &positions[0,0],
                   &mapping[0],
                   &n_classes,
                   &n_types,
                   &_symbols[0,0],
                   &atomic_numbers[0],
                   &occupations[0],
                   &eq_classes[0],
                   &n_symops,
                   &point_symmetry[0,0] if n_symops else NULL,
                   &point_symmetry[1,0] if n_symops else NULL,
                   &_mesh[0],
                   &n_symmetry_ops,
                   &rotations[0,0,0] if n_symmetry_ops >= 0 else NULL,
                   &translations[0,0] if n_symmetry_ops >= 0 else NULL,
                   &_verbose
                  )
    ratio = 1 / ratio
    centres[:max_spheres] *= ratio
    radii[:max_spheres + len(atoms)] *= Bohr
    spheres = Atoms(cell=atoms.cell, pbc=atoms.pbc, positions=centres[:max_spheres])
    spheres.wrap()
    out = EmptySpheresResult(spheres.positions, radii[:max_spheres])
    if return_atom_rws:
      return out, radii[max_spheres:max_spheres+len(atoms)]
    else:
      return out
