""" More atomic types can be on one site. """
from __future__ import annotations

import numpy as np
from .atomic_types import AtomicType
from typing import Dict, Union, Optional
from collections.abc import Iterable
from .sites import Site


class Occupation:
  """ Occupation of the atomic site, given by AtomicType : value

   The value determine the probability, that a given atomic type will
   be found on a given place (which can be used e.g. for computing
   alloys).
  """

  def __init__(self, dct:Dict[Union[AtomicType,str], float],
                     site:Optional[Site]=None, update_atoms=False):
      self._site = site
      self.set(dct, update_atoms)

  def copy(self, site:Optional[Site]=None) -> Occupation:
      """ Create a copy of the object, associated with a given site. """
      return Occupation({ a.copy(): v for a,v in self._occupation.items() }, site)

  def set(self, dct: Dict[AtomicType | str, float ], update_atoms=True) -> None:
      """ Set (replace) the occupation data.

      The method automatically updates the symbols, atomic numbers and the occupancy property (if exists)
      of the underlying Atoms object.
      """
      if isinstance(dct, (int, str, AtomicType)):
         dct = {dct : 1.0}
      if hasattr(dct, 'items'):
         iterator = dct.items()  # dict
      else:
         iterator = dct
      self._occupation = dict( (AtomicType.to_atomic_type(i), j) for i,j in iterator )
      self._normalize()
      if update_atoms:
          self._update_atoms()

  def items(self):
      """ dict.items() like enumeration """
      return self._occupation.items()

  def __repr__(self):
      return f"Occupation {self._occupation}"

  def __str__(self):
      return f"Occupation {self._occupation}"

  def _update_atoms(self):
      if self._site:
         self._site.update_atoms()

  def __iter__(self):
      return iter(self._occupation)

  def _find_key(self, name):
      if isinstance(name, int):
          return list(self._occupation.keys())[name]
      if isinstance(name, AtomicType):
          return name
      for i in self:
          if i.symbol == name:
              return i
      raise KeyError(f"No {name} in the occupation")

  def atomic_type(self, name: Union[str,int,AtomicType]) -> AtomicType | None:
      """ Find the corresponding atomic type according to the provided argument.

      Parameters
      ----------
      name
          The identification of the atomic type.
          If it is integer, returns the n-th atomic type (according to the orderd supplied
          when the occupation is set).
          If it is a string, the first atomic type of given chemical element is returned
          If it is AtomicType, it is returned "as is"
      """
      return self._find_key(name)

  def __getitem__(self, name):
      name = self._find_key(name)
      return self._occupation[name]

  def replace_type(self, name:Union[str,int,AtomicType], to:Union[str,AtomicType]):
      """
      Replace the given atomic type (see :meth:`atomic_type<ase2sprkkr.sprkkr.occupations.Occupation.atomic_type>`, how
      it can be identified) by the new one (given either by AtomicType or by its chemical symbol)
      """
      key = self._find_key(name)
      to = AtomicType.to_atomic_type(to)
      self._occupation = dict(
          (k if k is not key else to, v) for k,v in self._occupation.items()
        )

  def clean(self):
      """ Remove all items with zero probability. """
      self._occupation = dict(
          (k, v) for k,v in self._occupation.items() if v > 0
      )

  def __setitem__(self, name, value):
       try:
           name = self._find_key(name)
           self._normalize(1. - value, name)
           self._occupation[name] = value
       except KeyError:
           if isinstance(name, int):
              raise ValueError(f"Cannot add an atomic type using an integer key {name}. Please use a string (chemical symbol).")
           self.add(name, value)
       self._update_atoms()

  @property
  def primary_atomic_type(self):
      """ Return the atomic type.
          If there are more atoms on the site, return the one
          with the largest occupation """
      m = 0.
      prim = None
      for at,occ in self._occupation.items():
          if occ > m and at.atomic_number > 0:
             m = occ
             prim = at
      return prim

  @property
  def primary_atomic_number(self):
      """ Return the atomic number of the atom at the site.
          If there are more atoms on the site, return the "main one". """
      primary = self.primary_atomic_type
      return primary.atomic_number if primary else 0

  @property
  def primary_symbol(self):
      """ Return the chemical symbol of the atom at the site.
          If there are more atoms on the site, return the "main one". """
      primary = self.primary_atomic_type
      return primary.symbol if primary else 'X'

  def __len__(self):
      return len(self._occupation)

  def add(self, name, value=None):
      """ Add atom to the site """
      if value is None:
         value = 1. / (len(self) + 1.)
      self._normalize(1. - value)
      self._occupation[AtomicType.to_atomic_type(name)] = value
      self._update_atoms()

  def __delitem__(self, name):
      name = self._find_key(name)
      del self._occupation[name]
      self._normalize()

  def _normalize(self, to=1., except_from=None):
      """
      Normalizes occupation so the sum will be equal to the value 'to' (by default to 1.).
      If there are no None values, all the values are multiplied by the same number to make
      their sum equal to to.
      If there are None values, the remainder to the 'to' value is equally divided among them.

      Parameters
      ----------
      value : float
        Desired value

      except_from: AtomicType
        Skip the given atom during normalizing
      """
      suma = 0.
      none = []
      for i in self._occupation:
          if i is except_from:
             continue
          v = self._occupation[i]
          if v is None:
             none.append(i)
          else:
            suma+=v

      if suma >= to and none:
         for i in none:
             self._occupation[i] = 0.
         none = []
      if suma == to:
         return
      if none:
         for i in none:
             self._occupation[i] = (to - suma) / len(none)
      elif suma > to or except_from:
         ratio = to / suma

         for i in self._occupation:
             if i != except_from:
                 self._occupation[i] *= ratio
      else:
         for i in self._occupation:
            if i.is_vacuum():
               self._occupation[i]+= to - suma
         else:
               self._occupation[AtomicType('Vc')] = to - suma

  @property
  def as_dict(self):
      occ = {}
      for at in self:
         if at.atomic_number == 0:
            continue
         occ[at.symbol] = occ.get(at.symbol, 0) + self[at]
      return occ

  @as_dict.setter
  def as_dict(self, x):
      self.set(x)

  @property
  def total_occupation(self):
      return sum(self._occupation.values())

  @staticmethod
  def to_occupation(occupation, site):
      """ Create an occupation object associated with the given sites object """
      if not isinstance(occupation, Occupation):
         occupation = Occupation(occupation, site)
      elif site is not None:
         occupation = Occupation.copy(occupation, site)
      return occupation

  def check(self):
      if not np.isclose(sum(self),1.0):
         raise ValueError("Total occupation of the site should be equal to one")

  def to_tuple(self):
      return zip(self.keys(), self.values())

  def atomic_types(self) -> Iterable[AtomicType]:
      """ Returns the atomic types that can be present on the site. """
      return self._occupation.keys()
