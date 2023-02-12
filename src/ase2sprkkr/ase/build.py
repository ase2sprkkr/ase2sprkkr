import numpy as np
import ase
from ase.build import stack as ase_stack

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
          axis:int=None,
          at:Optional[List[Optional[List[int]]]]=None,
          relative:bool=False):
    """
    Stack (concatenate) the atoms objects along given axis

    Parameters
    ----------
    atomses
      List of atoms objects to be concatenated.

    axis
      Along which axis should be the atoms concatenated.
      The atoms are then stacked so the [0,0,0] relative coordinates
      of the (n+1)th atoms are located at [1,0,0], [0,1,0] or [0,0,1]
      respepectively (according to the axis argument) relative coordinates
      of the nth atoms object.
      If axis is not given, the at have to be supplied.

    at
      Determines the positions of the origins of coordinates of the atoms
      objects in the resulting objects. If it is None or [0,0,0] for the
      (n+1)th atoms, then the coordinates are determined as coordinates
      of the nth atoms plus axis-th cell vector.
      There can be n+1 items in the stack, then the last one determine the
      axis-th cell vector of the resulting atoms object.

    relative
      If True, the coordinates in at are considered as relative (to the coordinates of
      the previous atoms object).
    """
    out = atomses[0].copy()
    if at is None:
       valid_at = lambda n: False
    else:
       atlen = len(at)
       valid_at = lambda n: n<atlen and at[n] is not None and not np.equal(at[n],[0,0,0]).all()

    ax = axis or 0
    if valid_at(0):
       out.positions+=at[0]
       origin = at[0]
    else:
       origin = np.array([0.,0.,0.])

    shift = out.cell[ax]

    def update_origin(i):
       nonlocal origin
       if valid_at(i):
          if relative:
              origin+=at[i]
          else:
              origin=at[i]
       else:
          origin+=shift

    latoms = len(atomses)
    for i,a in zip(range(1,latoms), atomses[1:]):
       a = a.copy()
       out.pbc *= a.pbc
       update_origin(i)
       a.positions += origin
       out+=a
       shift=a.cell[ax]

    if axis:
      update_origin(latoms)
      out.cell[axis] = origin

    return out

def stack(atoms, atoms2=None, *args, **kwargs):

    if atoms2 is not None:
       return stack(atoms, atoms2, *args, **kwargs)
    else:
       out = atoms[0]
       for i in atoms[1:]:
           out = ase_stack(out, i, *args, **kwargs)
       return out

