""" This module contains routines for building materials.
Unlike ``ase2sprkkr.sprkkr.build``, this module contains generic
routines, possible usable with plain ASE (with any calculator).
"""

import numpy as np
import ase
from ase.build import stack as _stack

from typing import List, Union, Optional
def aperiodic_times(atoms:ase.Atoms,
                    times:Union[int, float, List[Union[int,float]]],
                    axis:Optional[int]=None,
                    direction:Union[List[int], int]=[1,1,1]):
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
    if isinstance(times, (int, float)):
       if axis is None:
          times = np.ones(3) * times
       else:
          t = np.ones(3)
          t[axis] = times
          times = t
    else:
       if axis is not None:
          raise ValueError("If axis is specified, only a scalar value have to "
                           "be supplied to the times argument.")

    if isinstance(direction, int):
       direction = [direction, direction, direction]

    assert len(times) == 3

    for i,(num, direction) in enumerate(zip(times, direction)):
        if times[i] == 1:
           continue
        inum = int(num)
        mlt = np.ones(3, dtype=int)
        mlt[i] = inum
        natoms = atoms*mlt
        fnum=num-inum
        if fnum>0:
           pos=atoms.get_scaled_positions(False)
           if direction > 0:
               add = atoms[pos[:,i] < fnum].copy()
               add.positions += atoms.cell[i] * inum
               natoms += add
           else:
               add = atoms[pos[:,i] > 1. - fnum].copy()
               add.positions -= atoms.cell[i]
               add += natoms
               natoms = add
               natoms.positions+= atoms.cell[i]*fnum
        natoms.cell[i] = atoms.cell[i]*num
        natoms.pbc[i]=False
        atoms = natoms
    return atoms

def stack(atomses:List[ase.Atoms],
          axis:int,
          at:Optional[List[Optional[List[int]]]]=None,
          relative:bool=False,
          scale='pbc',
          check_strain='auto',
          max_strain=1e-10,
          check_pbc=True,
          periodic=False):
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

    #first, define a function to retrieve the shifts
    if at is None:
       valid_at = lambda n: False
    else:
       atlen = len(at)
       valid_at = lambda n: n<atlen and at[n] is not None and not np.equal(at[n],[0,0,0]).all()

    def get_at(i):
        if not valid_at(i):
            return None
        a = at[i]
        if isinstance(a, (int, float)):
           out = out.cell[axis]
           out *= a / np.linalg.norm(out)
        else:
           out = a
        return out

    def update_origin(i):
       nonlocal origin
       a = get_at(i)
       if a is None:
          origin+=shift
       elif relative:
          origin+=a
       else:
          origin=a

    #set the initial origin and shift
    at0 = get_at(0)
    if at0 is None:
       origin = np.array([0.,0.,0.])
    else:
       out.positions+=at0
       origin = at0
    shift = out.cell[axis]


    #resolve resulting pbc
    cell_index = [ i for i in range(3) if i!=axis ]
    if check_strain == 'auto':
       check_strain = scale

    if not check_pbc:
        for a in remains:
            out.pbc *= a.pbc
    else:
        for a in remains:
            if (out.pbc != a.pbc)[cell_index].any():
                 raise ValueError("The stacked atoms has incompatibile pbc. Check the check_pbc argument.")
    out.pbc[axis] = periodic

    #and finally, stack the atoms
    a0cell = atoms0.cell.complete()
    for i,a in enumerate(remains, start=1):
       update_origin(i)
       out+=a
       out.pbc *= a.pbc
       positions = out.positions[-len(a):]

       #scaling of the incompatibile cells
       do_scale = []
       for c in cell_index:
          if (a.cell[c] != atoms0.cell[c]).any():
             if out.pbc[c] if check_strain == 'pbc' else check_strain:
                strain = np.linalg.norm(a.cell[c] - a0cell[c]) / np.linalg.norm(a0cell[c])
                if strain > max_strain:
                   raise ValueError("The {i}th stacked Atoms object {a.symbols} has incompatibile cell, check the max_strain argument.")
             if out.pbc[c] if scale == 'pbc' else scale:
                do_scale.append(c)
       if do_scale:
          cell = a.cell.complete()
          ncell = cell.copy()
          for c in do_scale:
              #copied from atoms.set_cell(scale_atoms=True)
              ncell[c] = a0cell[c]
          m = np.linalg.solve(cell, ncell)
          positions[:] = np.dot(positions, m)

       positions += origin
       shift=a.cell[axis]

    #update the cell of the resulting atoms
    update_origin(len(atomses))
    out.cell[axis] = origin
    return out
