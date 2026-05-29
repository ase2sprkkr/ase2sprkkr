"""This module contains routines for building materials.
Unlike ``ase2sprkkr.sprkkr.build``, this module contains generic
routines, possible usable with plain ASE (with any calculator).
"""

from fractions import Fraction
from math import gcd, lcm
import numbers

import numpy as np
import ase
from typing import List, Union, Optional, Tuple
from ase.build import surface


def aperiodic_times(
    atoms: ase.Atoms,
    times: Union[int, float, List[Union[int, float]]],
    axis: Optional[int] = None,
    direction: Union[List[int], int] = [1, 1, 1],
):
    """
    Multiply (repeat) the atoms in the same way as
    atoms __mult__ operator.

    However, it accepts floats, too, to add only a part of a last repeated cell.

    Parameters
    ----------
    atoms:
     The atoms to be repeated

    times:
     List of three integers or floats, that say, how many times should be
     the atoms repeated along given axis. If only one integer/float is given,
     that it is used for all axes, unless the axis argument is given.
     Float with a decimal part means, that the last (or first, see the direction
     argument) cell won't be added as whole, but only its part will be added.

    axis:
     If it is not None, the atoms are repeated

    direction:
     Integer or list of integers (one for each axis). If it is negative, the
     partial cell will be added on the begin of stacked cell, for the
     given axis. Otherwise, the partial cell

    """
    if isinstance(times, numbers.Real):
        if axis is None:
            times = np.ones(3) * times
        else:
            t = np.ones(3)
            t[axis] = times
            times = t
    else:
        if axis is not None:
            raise ValueError("If axis is specified, only a scalar value have to be supplied to the times argument.")

    if isinstance(direction, int):
        direction = [direction, direction, direction]

    assert len(times) == 3

    for i, (num, direction) in enumerate(zip(times, direction)):
        if times[i] == 1:
            continue
        inum = int(num)
        mlt = np.ones(3, dtype=int)
        mlt[i] = inum
        natoms = atoms * mlt
        fnum = num - inum
        if fnum > 0 and not np.isclose(fnum, 0):
            pos = atoms.get_scaled_positions(False)
            if direction > 0:
                add = atoms[pos[:, i] < fnum].copy()
                add.positions += atoms.cell[i] * inum
                natoms += add
            else:
                add = atoms[pos[:, i] > 1.0 - fnum].copy()
                add.positions -= atoms.cell[i]
                add += natoms
                natoms = add
                natoms.positions += atoms.cell[i] * fnum
            natoms.pbc[i] = False
        natoms.cell[i] = atoms.cell[i] * num
        atoms = natoms
    return atoms


def stack(
    atomses: List[ase.Atoms],
    axis: int,
    at: Optional[List[Optional[List[int]]]] = None,
    relative: bool = False,
    scale="pbc",
    check_strain="auto",
    max_strain=1e-10,
    check_pbc=True,
    periodic=False,
):
    """
    Stack (concatenate) the atoms objects along given axis

    This function is very similiar to ase.build.stack, but it
    support more than two atoms object to be stacked on themselves,
    and the arguments are a bit different.

    #TODO - could be this and ASE stack function merged together?


    Parameters
    ----------
    atomses
      List of atoms objects to be concatenated.

    axis
      Along which axis should be the atoms concatenated.
      The atoms are then stacked so the [0,0,0] relative coordinates
      of the (n+1)th atoms are located at [1,0,0], [0,1,0] or [0,0,1]
      respepectively (according to the axis argument) relative cell
      coordinates of the nth atoms object.
      The at can shift these distances.

    at
      Determines the positions of the origins of coordinates of the atoms
      objects in the resulting objects. If it is None or [0,0,0] for the
      (n+1)th atoms, then the coordinates are determined as the coordinates
      of the nth atoms plus the axis-th cell vector.
      There can be n+1 items in the stack, then the last one determine the
      axis-th cell vector of the resulting atoms object.
      If the given item is just one scalar r, it is considered as r*unitary
      vector along the axis.

    relative
      If True, the coordinates in at are considered as relative to the
      cell corner (see the axis argument).

    scale
      If True, the stacked atoms are scaled in the two dimensions (not
      in the axis one) so the corresponding two cell vectors
      are the same as these of the first atoms.

      The default value ``'pbc'`` means, that scaling is done only if the
      given axis is periodic.

    check_strain
      Check the compatibility of the cells along the other two
      axes (not the one along which the atoms are stacked).

      If True, the maximal strain cannot exceeded the max_strain argument.
      The default value ``'auto'`` means the same value as scale.
      The value ``'pbc'`` means check the strain only along the axes that
      are periodic.

      If the strain is exceeded, a ValueError is raised.

    max_strain
      The limit for a maximal (relative to the norm of the corresponding
      first atoms cell vector) displacement of the cell vectors.

    check_pbc
      If True, all the atoms objects have to have the same pbc along the other two axes.

    periodic
      The pbc of the resulting object along the axis.
    """
    try:
        atoms0 = atomses[0]
        remains = atomses[1:]
    except TypeError:
        iterator = iter(atomses)
        atoms0 = next(iterator)
        remains = [i for i in iterator]

    out = atoms0.copy()

    # first, define a function to retrieve the shifts
    if at is None:

        def valid_at(n):
            return False
    else:
        atlen = len(at)

        def valid_at(n):
            return n < atlen and at[n] is not None and not np.equal(at[n], [0, 0, 0]).all()

    def get_at(i):
        if not valid_at(i):
            return None
        a = at[i]
        if isinstance(a, (int, float)):
            o = out.cell[axis]
            o *= a / np.linalg.norm(out)
        else:
            o = a
        return o

    def update_origin(i):
        nonlocal origin
        a = get_at(i)
        if a is None:
            origin += shift
        elif relative:
            origin += a
        else:
            origin = a

    # set the initial origin and shift
    at0 = get_at(0)
    if at0 is None:
        origin = np.array([0.0, 0.0, 0.0])
    else:
        out.positions += at0
        origin = at0
    shift = out.cell[axis]

    # resolve resulting pbc
    cell_index = [i for i in range(3) if i != axis]
    if check_strain == "auto":
        check_strain = scale

    if not check_pbc:
        for a in remains:
            out.pbc *= a.pbc
    else:
        for a in remains:
            if (out.pbc != a.pbc)[cell_index].any():
                raise ValueError("The stacked atoms has incompatibile pbc. Check the check_pbc argument.")
    out.pbc[axis] = periodic

    # and finally, stack the atoms
    a0cell = atoms0.cell.complete()
    for i, a in enumerate(remains, start=1):
        update_origin(i)
        out += a
        out.pbc *= a.pbc
        positions = out.positions[-len(a) :]

        # scaling of the incompatibile cells
        do_scale = []
        for c in cell_index:
            if (a.cell[c] != atoms0.cell[c]).any():
                if out.pbc[c] if check_strain == "pbc" else check_strain:
                    strain = np.linalg.norm(a.cell[c] - a0cell[c]) / np.linalg.norm(a0cell[c])
                    if strain > max_strain:
                        raise ValueError(
                            "The {i}th stacked Atoms object {a.symbols} has incompatibile cell, check the max_strain argument."
                        )
                if out.pbc[c] if scale == "pbc" else scale:
                    do_scale.append(c)
        if do_scale:
            cell = a.cell.complete()
            ncell = cell.copy()
            for c in do_scale:
                # copied from atoms.set_cell(scale_atoms=True)
                ncell[c] = a0cell[c]
            m = np.linalg.solve(cell, ncell)
            positions[:] = np.dot(positions, m)

        positions += origin
        shift = a.cell[axis]

    # update the cell of the resulting atoms
    update_origin(len(atomses))
    out.cell[axis] = origin
    return out


def _surface_basis(cell, indices, tol=1e-10):
    indices = np.asarray(indices)
    if indices.shape != (3,) or not indices.any():
        raise ValueError(f"{indices} is an invalid surface type")
    if not np.allclose(indices, np.rint(indices), atol=tol):
        raise ValueError(f"{indices} is an invalid surface type")

    indices = np.rint(indices).astype(int)
    h, k, l = indices
    h0, k0, l0 = indices == 0

    if h0 and k0 or h0 and l0 or k0 and l0:
        if not h0:
            c1, c2, c3 = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
        if not k0:
            c1, c2, c3 = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
        if not l0:
            c1, c2, c3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    else:
        p, q = _ext_gcd(k, l)
        a1, a2, a3 = np.asarray(cell)

        k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3), l * a2 - k * a3)
        k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3), l * a2 - k * a3)

        if abs(k2) > tol:
            i = -int(round(k1 / k2))
            p, q = p + i * l, q - i * k

        a, b = _ext_gcd(p * k + q * l, h)

        c1 = (p * k + q * l, -p * h, -q * h)
        c2 = np.array((0, l, -k)) // abs(gcd(l, k))
        c3 = (b, a * p, a * q)

    return np.array([c1, c2, c3], dtype=float)


def _ext_gcd(a, b):
    if b == 0:
        return 1, 0
    if a % b == 0:
        return 0, 1
    x, y = _ext_gcd(b, a % b)
    return y, x - y * (a // b)


def _fractional_inplane_period(shift: np.ndarray, tol: float = 1e-10, max_denominator: int = 256) -> Optional[int]:
    period = 1
    for value in np.asarray(shift, dtype=float):
        fraction = value - np.floor(value)
        if np.isclose(fraction, 0.0, atol=tol) or np.isclose(fraction, 1.0, atol=tol):
            continue
        rational = Fraction(fraction).limit_denominator(max_denominator)
        if abs(float(rational) - fraction) > tol:
            return None
        period = lcm(period, rational.denominator)
    return period


def _minimal_surface_layers(atoms: ase.Atoms, hkl: Tuple[float], tol: float = 1e-10, max_layers: int = 256) -> int:
    basis = _surface_basis(atoms.cell, hkl, tol=tol)
    cell = np.dot(basis, np.asarray(atoms.cell))
    a1, a2, a3 = cell

    metric = np.array([[np.dot(a1, a1), np.dot(a1, a2)], [np.dot(a2, a1), np.dot(a2, a2)]])
    rhs = np.array([np.dot(a1, a3), np.dot(a2, a3)])
    shift = np.linalg.solve(metric, rhs)

    period = _fractional_inplane_period(shift, tol=tol, max_denominator=max_layers)
    if period is not None:
        return period

    for layers in range(1, max_layers + 1):
        if np.allclose(layers * shift, np.rint(layers * shift), atol=tol):
            return layers

    raise ValueError(f"Unable to determine a periodic slab thickness for Miller indices {tuple(hkl)}.")


def rotate(atoms: ase.Atoms, hkl: Tuple[float]):
    """
    Rotate the atoms according the given Miller coordinates.

    Parameters
    ----------
    atoms
      The atoms to be rotated and shifted
    hkl
      Miller indices. The atoms will be rotated so that the last axis will be perpendicular
      to the plane and the other will be parallel.
    """
    layers = _minimal_surface_layers(atoms, hkl)
    return surface(atoms, np.rint(hkl).astype(int), layers, periodic=True)


def shift(atoms: ase.Atoms, shift: Optional[Union[float, int, tuple, list, np.ndarray]], axis=2, wrap=True):
    """
    Shift the atoms (to get the desired atom to the top/bottom of the cell).

    Parameters
    ----------
    atoms
      The atoms to shift their positions

    shift
      If a vector is given, add it to the positions.
      If a float is given, shift byt the given fraction of the given axis
      If an integer is given, shift so that the given atom is at [0,0,0]
    """
    if isinstance(shift, int):
        shift = -atoms.positions[shift]
    else:
        shift = atoms.cell[axis] * shift
    atoms.positions[:] = atoms.positions + shift
    if wrap:
        atoms.wrap()
    return atoms


def flip_around(atoms, around, axis=2, wrap=True):
    """
    Flip the atoms in given axis around given point (or site)
    so that the point/site will be the most distant from origin of all the sites.
    """
    shift(atoms, -around, wrap=False)
    flip(atoms, axis, wrap=False)
    sp = atoms.get_scaled_positions()
    aid = np.argmax(sp[:, axis] % 1.0)
    return shift(atoms, sp[aid, axis], wrap=wrap)


def flip(atoms, plane=2, wrap=True):
    """
    Flip the atoms around plane perpendicular to the given vector (if vector is given)
    or axis (identified by an integer)
    """
    if isinstance(plane, int):
        vc = atoms.cell[plane]
    else:
        vc = plane
    vc = vc / np.linalg.norm(vc)
    reflect = np.eye(3) - 2 * np.outer(vc, vc)
    atoms.positions = np.dot(atoms.positions, reflect.T)

    if wrap:
        atoms.wrap()
    return atoms
