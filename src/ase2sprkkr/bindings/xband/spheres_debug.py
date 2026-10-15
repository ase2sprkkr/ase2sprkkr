from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from ase import Atoms
from ase.units import Bohr

from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ...common.unique_values import UniqueValuesMapping


@dataclass(frozen=True)
class FindEmptySpheresInputs:
    max_spheres: int
    min_radius_bohr: float
    max_radius_bohr: float
    alat_bohr: float
    cell_scaled: np.ndarray  # shape (3,3), C-order like in spheres.pyx
    positions_scaled: np.ndarray  # shape (n,3)
    mapping: np.ndarray  # shape (n,), int32
    n_classes: int
    n_types: int
    symbols4: list[str]  # len n_types, each exactly 4 chars
    atomic_numbers: np.ndarray  # shape (n_types,), float64
    occupations: np.ndarray  # shape (n_types,), float64
    type_eq_class: np.ndarray  # shape (n_types,), int32
    mesh: np.ndarray  # shape (3,), int32
    n_symmetry_ops: int
    rotations: np.ndarray | None  # shape (n_symmetry_ops, 3,3)
    translations: np.ndarray | None  # shape (n_symmetry_ops, 3)
    verbose: int


def _fmt_f64(x: float) -> str:
    # Force exponent form that gfortran accepts and is stable for copy-paste.
    return f"{float(x):.17e}".replace("e", "d")


def _fmt_i32(x: int) -> str:
    return str(int(x))


def _pad4(s: str) -> str:
    s = (s or "")[:4]
    return (s + "    ")[:4]


def _as_list_literal_f64(values: Iterable[float]) -> str:
    return ", ".join(_fmt_f64(v) for v in values)


def _as_list_literal_i32(values: Iterable[int]) -> str:
    return ", ".join(_fmt_i32(v) for v in values)


def _collect_inputs(
    atoms: Atoms,
    *,
    min_radius: float,
    max_radius: float,
    mesh: int | Sequence[int],
    verbose: bool | int,
    max_spheres: int,
) -> FindEmptySpheresInputs:
    SPRKKRAtoms.promote_ase_atoms(atoms)

    to_bohr = 1.0 / Bohr
    min_radius_bohr = float(min_radius) * to_bohr
    max_radius_bohr = float(max_radius) * to_bohr

    alat_bohr = atoms.cell.get_bravais_lattice().a * to_bohr

    es = atoms.spacegroup_info.equivalent_sites
    es = UniqueValuesMapping(es)
    ui = es.unique_indexes()
    mapping = np.asarray(es.normalized(dtype=np.int32)[0], dtype=np.int32)

    n_types = sum(len(atoms.sites[i].occupation) for i in ui)
    n_classes = len(ui)

    symbols4: list[str] = []
    occupations = np.empty(n_types, dtype=np.float64)
    atomic_numbers = np.empty(n_types, dtype=np.float64)
    type_eq_class = np.empty(n_types, dtype=np.int32)

    type_no = 0
    for index in ui:
        site = atoms.sites[index]
        for typ, occ in site.occupation.items():
            symbols4.append(_pad4(getattr(typ, "symbol", str(typ))))
            occupations[type_no] = float(occ)
            atomic_numbers[type_no] = float(getattr(typ, "atomic_number", 0.0))
            type_eq_class[type_no] = int(index) + 1  # 1-based, matches spheres.pyx
            type_no += 1

    if isinstance(mesh, int):
        mesh3 = np.array([mesh, mesh, mesh], dtype=np.int32)
    else:
        mesh3 = np.asarray(mesh, dtype=np.int32)
        if mesh3.shape != (3,):
            raise ValueError("mesh must be int or length-3 sequence")

    if isinstance(verbose, bool):
        verbose_i = 1 if verbose else 0
    else:
        verbose_i = int(verbose)

    ratio = to_bohr / float(alat_bohr)
    cell_scaled = np.asarray(atoms.cell[:], dtype=np.float64) * ratio
    positions_scaled = np.asarray(atoms.positions, dtype=np.float64) * ratio

    sp = atoms.spacegroup_info
    rotations = None
    translations = None
    n_symmetry_ops = -1

    if sp and sp.dataset and sp.dataset.rotations is not None:
        d = sp.dataset
        rotations = np.asarray(d.rotations, dtype=np.float64)  # shape (nops,3,3)
        translations = np.asarray(d.translations, dtype=np.float64)  # shape (nops,3)
        n_symmetry_ops = int(len(rotations))

        # Emit rotations exactly as spglib provides them (nops,3,3) but stored
        # into a Fortran-declared array ROTATIONS(3,3,nops) in the same numeric
        # orientation: ROTATIONS(i,j,op) = rot_py(op,i,j).

    return FindEmptySpheresInputs(
        max_spheres=int(max_spheres),
        min_radius_bohr=min_radius_bohr,
        max_radius_bohr=max_radius_bohr,
        alat_bohr=float(alat_bohr),
        cell_scaled=cell_scaled,
        positions_scaled=positions_scaled,
        mapping=mapping,
        n_classes=int(n_classes),
        n_types=int(n_types),
        symbols4=symbols4,
        atomic_numbers=atomic_numbers,
        occupations=occupations,
        type_eq_class=type_eq_class,
        mesh=mesh3,
        n_symmetry_ops=int(n_symmetry_ops),
        rotations=rotations,
        translations=translations,
        verbose=int(verbose_i),
    )


def fortran_code_find_empty_spheres(
    atoms: Atoms,
    *,
    min_radius: float = 0.65,
    max_radius: float = 2.0,
    mesh: int | Sequence[int] = 24,
    verbose: bool | int = False,
    max_spheres: int = 256,
    program_name: str = "repro_find_empty_spheres",
) -> str:
    """Return a standalone Fortran program that calls `FIND_EMPTY_SPHERES`.

    The emitted program mirrors `spheres.pyx` data preparation:
    - radii in Bohr
    - `alat` in Bohr
    - `CELL` and `BAS` scaled by `ratio = (1/Bohr)/alat`
    - symmetry rotations/translations layout as seen by Fortran

    It is meant for debugging crashes / wrong results by reproducing the exact call.
    """

    inp = _collect_inputs(
        atoms, min_radius=min_radius, max_radius=max_radius, mesh=mesh, verbose=verbose, max_spheres=max_spheres
    )

    n = int(inp.positions_scaled.shape[0])
    n_out = inp.max_spheres

    # Fortran expects NTMAX=200 in radii.f90; keep that as a conservative max.
    ntmax = 200
    if inp.n_types > ntmax:
        raise ValueError(f"n_types={inp.n_types} exceeds NTMAX={ntmax} used by FIND_EMPTY_SPHERES")
    if n > ntmax:
        raise ValueError(f"n_atoms={n} exceeds NTMAX={ntmax} used by FIND_EMPTY_SPHERES")

    lines: list[str] = []
    a = lines.append

    a(f"program {program_name}")
    a("  implicit none")
    a("  integer, parameter :: NTMAX=200")
    a("  integer, parameter :: MAXDIM=1000000")
    a("")
    a("  integer :: ret")
    a("  integer :: n_out")
    a("  integer :: nq")
    a("  integer :: nm")
    a("  integer :: nsort")
    a("  integer :: verbose")
    a("  integer :: n_symmetry_ops")
    a("  real(8) :: rmines, rmaxes, alat")
    a("  integer :: mesh(3)")
    a("  real(8) :: cell(3,3)")
    a("  real(8) :: centres(3, MAXDIM)")
    a("  real(8) :: radii(MAXDIM)")
    a("  real(8) :: bas(3, NTMAX)")
    a("  integer :: imq(NTMAX)")
    a("  integer :: imt(NTMAX)")
    a("  character(len=4) :: txtt(NTMAX)")
    a("  real(8) :: z(NTMAX)")
    a("  real(8) :: conc(NTMAX)")
    a("  real(8) :: rotations(3,3,96)")
    a("  real(8) :: translations(3,96)")
    a("  integer :: i")
    a("")

    a("  integer, external :: find_empty_spheres")
    a("")

    a(f"  n_out = {n_out}")
    a(f"  nq = {n}")
    a(f"  nm = {inp.n_classes}")
    a(f"  nsort = {inp.n_types}")
    a(f"  verbose = {inp.verbose}")
    a(f"  rmines = {_fmt_f64(inp.min_radius_bohr)}")
    a(f"  rmaxes = {_fmt_f64(inp.max_radius_bohr)}")
    a(f"  alat = {_fmt_f64(inp.alat_bohr)}")
    a(f"  mesh = (/ {_as_list_literal_i32(inp.mesh.tolist())} /)")

    a("  cell(:,:) = 0d0")
    # Cell rows -> Fortran columns (matches spheres.pyx + Fortran column-major view)
    for col in range(3):
        row = inp.cell_scaled[col, :]
        a(
            "  cell(1:3,{j}) = (/ {x0}, {x1}, {x2} /)".format(
                j=col + 1, x0=_fmt_f64(row[0]), x1=_fmt_f64(row[1]), x2=_fmt_f64(row[2])
            )
        )

    a("  bas(:,:) = 0d0")
    for i in range(n):
        p = inp.positions_scaled[i, :]
        a(
            "  bas(1:3,{j}) = (/ {x0}, {x1}, {x2} /)".format(
                j=i + 1, x0=_fmt_f64(p[0]), x1=_fmt_f64(p[1]), x2=_fmt_f64(p[2])
            )
        )

    a("  imq(:) = 0")
    # mapping is 1-based already
    for i in range(n):
        a(f"  imq({i + 1}) = {_fmt_i32(int(inp.mapping[i]))}")

    a("  txtt(:) = '    '")
    a("  z(:) = 0d0")
    a("  conc(:) = 0d0")
    a("  imt(:) = 0")
    for i in range(inp.n_types):
        a(f"  txtt({i + 1}) = '{inp.symbols4[i]}'")
        a(f"  z({i + 1}) = {_fmt_f64(float(inp.atomic_numbers[i]))}")
        a(f"  conc({i + 1}) = {_fmt_f64(float(inp.occupations[i]))}")
        a(f"  imt({i + 1}) = {_fmt_i32(int(inp.type_eq_class[i]))}")

    a("  rotations(:,:,:) = 0d0")
    a("  translations(:,:) = 0d0")

    if inp.n_symmetry_ops >= 0 and inp.rotations is not None and inp.translations is not None:
        if inp.n_symmetry_ops > 96:
            raise ValueError(f"n_symmetry_ops={inp.n_symmetry_ops} exceeds 96 hardcoded in emitted program")
        a(f"  n_symmetry_ops = {inp.n_symmetry_ops}")
        for op in range(inp.n_symmetry_ops):
            r = inp.rotations[op]
            t = inp.translations[op]
            op += 1
            # rotations(:, :, op+1)
            for j in range(3):
                a(
                    "  rotations(:,{j}, {op}) = (/ {x0}, {x1}, {x2} /)".format(
                        op=op, j=j + 1, x0=_fmt_f64(r[j, 0]), x1=_fmt_f64(r[j, 1]), x2=_fmt_f64(r[j, 2])
                    )
                )
            a(
                "  translations(1:3,{op}) = (/ {x0}, {x1}, {x2} /)".format(
                    op=op, x0=_fmt_f64(t[0]), x1=_fmt_f64(t[1]), x2=_fmt_f64(t[2])
                )
            )
    else:
        a("  n_symmetry_ops = -1")

    a("")
    a("  ret = find_empty_spheres( &")
    a("    n_out, centres, radii, rmines, rmaxes, alat, cell, nq, bas, &")
    a("    imq, nm, nsort, txtt, z, conc, imt, mesh, &")
    a("    n_symmetry_ops, rotations, translations, verbose &")
    a("  )")
    a("  print *, 'FIND_EMPTY_SPHERES ret=', ret")
    a("  print *, 'N_OUT=', n_out")
    a("  DO I=1,N_OUT")
    a("     print *, CENTRES(:,I), RADII(I)")
    a("  END DO")
    a("end program")

    return "\n".join(lines) + "\n"


def print_find_empty_spheres_program(*args, **kwargs) -> None:
    print(fortran_code_find_empty_spheres(*args, **kwargs), end="")
